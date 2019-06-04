#!/usr/bin/env bash
#
# Processes the dta file for rRNA contamination


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

find ${R_RNA_PATH} -name "*.fasta" -maxdepth 1 |
    while read f;
    do
        name=$(echo ${f}| xargs -I {} basename {} ) # Get file name
        # Check if indexing has done already
        if [[ -n $(find "${R_RNA_PATH}/index/" -name "${name}*") ]]; then
            log "Indexing has already done for ${name}"
        else
            ${TOOL_SORTMERNA}/bin/indexdb --ref ${f},${R_RNA_PATH}/index/${name}.idx -v # Make Index
        fi
    done

log "rRNA indexing completed successfully and stored in '${R_RNA_PATH}/index/'"
log "Checking if merged file required for filtering exists"

readonly MERGE_FOLDER="${DATA}/filtering/merged"

mkdir -p ${MERGE_FOLDER}

# This analysis assume Paired-end SRA file which has split into 2 parts

readonly FORWARD_READS="${DATA}/${SRA_ID}_1.fastq"
readonly REVERSE_READS="${DATA}/${SRA_ID}_2.fastq"

if [[ -n $(find ${MERGE_FOLDER} -name "${SRA_ID}_merged.fastq") ]]; then
    log "Merged file already exists in '${MERGE_FOLDER}'"
else
    ${TOOL_SORTMERNA}/scripts/merge-paired-reads.sh ${FORWARD_READS} ${REVERSE_READS} "${MERGE_FOLDER}/${SRA_ID}_merged.fastq"
    log "Merged reads generated for rRNA filtering"
fi

