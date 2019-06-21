#!/usr/bin/env bash
#
# Processes the data file for rRNA contamination

source ./config.sh

log "Checking for rRNA contamination"

# Path of folder where you want your fasta files and index files to be stored
readonly R_RNA_PATH="$GENOME_HOME/rrna"

# Index files will be stored here
mkdir -p ${R_RNA_PATH}/index

# Check if rRNA FASTA files exists

readonly TEMP2=$(find ${R_RNA_PATH} -name "*.fasta"| wc -l)

# This will skip copying files from SortMeRNA even if there is single rRNA FASTA file present in this folder

if [[ ${TEMP2} -gt 0 ]]; then
    log "Total of ${TEMP2} rRNA FASTA files found in '${R_RNA_PATH}'"
    log "Skipping copying FATSA files from SortMeRNA database"
else
    log "No rRNA FATSA files found. Trying to copy from SortMeRNA database"

    # Assuming all rRNA database is in 'rRNA_databases' folder
    find "${TOOL_SORTMERNA}/rRNA_databases/" -name "*.fasta"|
    while read f;
    do
        cp ${f} ${R_RNA_PATH} # Copy to our reference folder
    done

    log "rRNA FATSA files copied from SortMeRNA database"
fi

log "Starting rRNA indexing"

find ${R_RNA_PATH} -maxdepth 1 -name "*.fasta" |
    while read f;
    do
        name=$(echo ${f}| xargs -I {} basename {} ) # Get file name
        # Check if indexing has done already
        if [[ -n $(find "${R_RNA_PATH}/index/" -name "${name}*") ]]; then
            log "Indexing has already done for ${name}"
        else
            ${TOOL_SORTMERNA}/bin/indexdb --ref ${f},${R_RNA_PATH}/index/${name}${R_RNA_INDEX_EXTENSION} -v # Make Index
        fi
    done

log "rRNA indexing completed successfully and stored in '${R_RNA_PATH}/index/'"
log "Checking if merged file required for filtering exists"

readonly MERGE_FOLDER="${DATA}/filtering/merged"

mkdir -p ${MERGE_FOLDER}

# This analysis assume Paired-end SRA file which has split into 2 parts

readonly FORWARD_READS="${DATA}/${SRA_ID}_1.fastq"
readonly REVERSE_READS="${DATA}/${SRA_ID}_2.fastq"

if [[ -n $(find ${MERGE_FOLDER} -name "${SRA_ID}${READ_MERGED_EXTENSION}") ]]; then
    log "Merged file already exists in '${MERGE_FOLDER}'"
else
    ${TOOL_SORTMERNA}/scripts/merge-paired-reads.sh \
    ${FORWARD_READS} \
    ${REVERSE_READS} \
    "${MERGE_FOLDER}/${SRA_ID}${READ_MERGED_EXTENSION}"

    log "Merged reads generated for rRNA filtering"
fi

# Generate sequence that can be used with SortMeRNA command

REF_SEQ=""
while read f;
do
    name=$(echo ${f}| xargs -I {} basename {} ) # Get file name
    REF_SEQ="${REF_SEQ}${f},${R_RNA_PATH}/index/${name}${R_RNA_INDEX_EXTENSION}:"
done <<< $(find ${R_RNA_PATH} -maxdepth 1 -name "*.fasta" -type f)

REF_SEQ=$(echo ${REF_SEQ} | sed 's/.$//g')

log "rRNA filtering initiating"

# Make folder to keep aligned sequences

FILTERED_FOLDER="${DATA}/filtering"
mkdir -p ${FILTERED_FOLDER}

# Make Folder for Key-value data-store

KVDS_PATH="${DATA}/filtering/kvds"
rm -rf ${KVDS_PATH} # This is needed because sortmerna script won't run if this folder is non-empty
mkdir -p ${KVDS_PATH}


# Start SortMeRNA
#
# --ref : rRNA reference sequences and their indices
# --reads : Merged Fasta file
# --aligned : Aligned file name and path
# --other : non-aligned file name and path (this is non-rRNA)
# --fastx : fasta format
# --log : Creates log files
# -v : Verbose output
# -a : number of threads (if you don't know how many core your system has, you can ommit this option)
# -d : key value data store dictionary path (This option is needed because sortmerna won't delete kvds folder file by itself)
# --paired_in : If either one from forward or reverse read matches, it will be in aligned file

# Check if Filtered sequence is present

if [[ -n $(find "${FILTERED_FOLDER}" -maxdepth 1 -type f -name "${SRA_ID}${R_RNA_FILTERED_EXTENSION}.fastq") ]]; then
    log "Filtered file already exist. Skipping rRNA filtering"
else
    log "starting rRNA filtering"

    ${TOOL_SORTMERNA}/bin/sortmerna \
    --ref ${REF_SEQ} \
    --reads "${MERGE_FOLDER}/${SRA_ID}${READ_MERGED_EXTENSION}" \
    --aligned "${FILTERED_FOLDER}/${SRA_ID}_rRNA" \
    --other "${FILTERED_FOLDER}/${SRA_ID}${R_RNA_FILTERED_EXTENSION}" \
    --fastx \
    --paired_in \
    --log \
    -v \
    -a 8 \
    -d ${KVDS_PATH}

    log "rRNA filtering completed"
fi

# Split the files



TEMP3=$(find ${FILTERED_FOLDER} -maxdepth 1 -type f -name "${SRA_ID}${R_RNA_FILTERED_EXTENSION}_*.fastq" | wc -l)

if [[ ${TEMP3} -eq 2 ]]; then
    log "${TEMP3} filtered-split sequences already exists in ${FILTERED_FOLDER}"
    log "Skipping splitting of filtered sequence"
else
    log "Splitting filtered sequence"

    ${TOOL_SORTMERNA}/scripts/unmerge-paired-reads.sh \
    "${FILTERED_FOLDER}/${SRA_ID}${R_RNA_FILTERED_EXTENSION}.fastq" \
    "${FILTERED_FOLDER}/${SRA_ID}${R_RNA_FILTERED_EXTENSION}_1.fastq" \
    "${FILTERED_FOLDER}/${SRA_ID}${R_RNA_FILTERED_EXTENSION}_2.fastq"

    log "Splitting completed in folder ${FILTERED_FOLDER}"
fi

log "Checking quality of filtered files"

# Make folder to keep analysis files
QUALITY_FOLDER="${DATA}/quality/${SRA_ID}/"
mkdir -p ${QUALITY_FOLDER}

# Get all the fastq files for analysis

# Check if quality analysis is already done
find ${FILTERED_FOLDER} -maxdepth 1 -name "${SRA_ID}${R_RNA_FILTERED_EXTENSION}_*.fastq" |
while read f;
do
    name=$(echo ${f}| xargs -I {} basename {} )
    check_name=$(echo ${name} | sed -r 's/\.fastq/\_/')
    if [[ -n $(find ${QUALITY_FOLDER} -name "${check_name}*") ]]; then
        log "Quality score already exists for ${name}. Skipping calculation of quality"
    else
        ${TOOL_FASTQC}/fastqc \
        ${f} \
        --outdir=${QUALITY_FOLDER} # Output Folder for generated result files
    fi
done

log "Quality checked for filtered FASTQ files. results can be found in '${QUALITY_FOLDER}'"
log "All files are ready for final alignment"
