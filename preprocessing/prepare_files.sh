#!/usr/bin/env bash
#
# This scripts downloads SRA file from NCBI server, converts it to fastq format and generates their quality score.
# If any of already analysed or downloaded files are present, analysis will automatically skip the step

source ./config.sh

log "Preparing Files"
log "Checking if file is already downloaded"

# Check if file is in raw folder
if [[ -f "${DATA}/raw/${SRA_ID}.sra" ]]; then
    log "File '${SRA_ID}.sra' already exists in raw data. Skipping download."
# Check file is in NCBI's local downloading folder
elif [[ -f "${NCBI_DOWNLOAD}/${SRA_ID}.sra" ]]; then
    log "${SRA_ID}.sra file found in ncbs local folder. Moving to raw folder"
    mv "${NCBI_DOWNLOAD}/${SRA_ID}.sra" ${DATA}/raw/
# If not, download with prefetch
else
    log "File '${SRA_ID}.sra' not found. Downloading with prefetch"
    ${TOOL_SRA}/bin/prefetch -c ${SRA_ID}
fi

log "File '${SRA_ID}.sra' is ready for further processing"

# Following command is borrowed from  https://edwards.sdsu.edu/research/fastq-dump/
#
#--gzip : Compression format (not required here)
#--skip-technical : Exclude technical reads and includes only biological reads
#--readids : Append IDs if fies are split by spot
#--dumpbase : Formats sequence using base space (only required for SOLiD)
#--split-files : Dump each read into separate file.
#--clip : Removes left and right clips/tags
#--outdir : If we want output in different dictionary
#PATH/TO/SRA_FILE : Full path of file

readonly TEMP1=$(find ${DATA} -name "${SRA_ID}*.fastq" -maxdepth 1 | wc -l)

log "Checking if FASTQ file exists"
if [[ ${TEMP1} -gt 0 ]]; then
    log "Total of ${TEMP1} FASTQ files exists. Skipping FASTQ conversion."
else
    log "FASTQ files not found. Starting SRA to FASTQ conversion"

    ${TOOL_SRA}/bin/fastq-dump \
    --skip-technical \
    --readids \
    --dumpbase \
    --split-files \
    --clip \
    --outdir ${DATA} \
    "${DATA}/raw/${SRA_ID}.sra"

    log "SRA files converted to FASTQ for further analysis"
fi

# Check Quality of the data
log "Checking generated FASTQ files for their quality"

# Make folder to keep analysis files
QUALITY_FOLDER="${DATA}/quality/${SRA_ID}/"
mkdir -p ${QUALITY_FOLDER}

# Get all the fastq files for analysis

# Check if quality analysis is already done
find ${DATA} -maxdepth 1 -name "${SRA_ID}*.fastq" |
while read f;
do
    name=$(echo ${f}| xargs -I {} basename {} )
    check_name=$(echo ${name} | sed "s/\./\_/g")
    if [[ -n $(find ${QUALITY_FOLDER} -name "${check_name}*") ]]; then
        log "Quality score already exists for ${name}. Skipping calculation of quality"
    else
        ${TOOL_FASTQC}/fastqc \
        ${f} \
        --outdir=${QUALITY_FOLDER} # Output Folder for generated result files
    fi
done

log "Quality checked for FASTQ files. results can be found in '${QUALITY_FOLDER}'"



