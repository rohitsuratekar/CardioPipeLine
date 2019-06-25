#!/usr/bin/env bash
#
# Alignment and Mapping using STAR
#
# --runThreadN : Number of threads to run
# --genomeDir : Path of folder where indices are generated. (In this pipeline from "star_indexing.sh" script)
# --readFilesIn : Reads to be mapped
# --outFileNamePrefix : Output file path

source ./config.sh

STAR_INDEX_PATH="index/star"
READ_1="${DATA}/filtering/${SRA_ID}${R_RNA_FILTERED_EXTENSION}_1.fastq"
READ_2="${DATA}/filtering/${SRA_ID}${R_RNA_FILTERED_EXTENSION}_2.fastq"
OUTPUT_PATH="${DATA}/aligned/${SRA_ID}"
mkdir -p ${OUTPUT_PATH}
OUTPUT_FILE="${OUTPUT_PATH}/${SRA_ID}_"

if [[ ${SKIP_RNA_FILTERING} -eq 1 ]]; then
    log "rRNA filtering is not done. Using file without rRNA filtering"
    READ_1="${DATA}/${SRA_ID}_1.fastq"
    READ_2="${DATA}/${SRA_ID}_2.fastq"
else
    log "Initializing STAR Alignment"
fi



${TOOL_STAR}/STAR \
--runThreadN 8 \
--genomeDir ${GENOME_HOME}/${STAR_INDEX_PATH} \
--readFilesIn ${READ_1} ${READ_2} \
--outFileNamePrefix ${OUTPUT_FILE}


log "STAR Alignment finished"