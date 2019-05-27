#!/usr/bin/env bash
#
# Alignment and Mapping using STAR
#
# --runThreadN : Number of threads to run
# --genomeDir : Path of folder where indices are generated. (In this pipeline from "star_indexing.sh" script)
# --readFilesIn : Reads to be mapped
# --outFileNamePrefix : Output file path

source ./config.sh

readonly STAR_INDEX_PATH="index/star"
readonly READ_1="${SRA_ID}_pass_1.fastq"
readonly READ_2="${SRA_ID}_pass_2.fastq"
readonly OUTPUT_PATH="${DATA}/aligned"
mkdir -p ${OUTPUT_PATH}
readonly OUTPUT_FILE="${OUTPUT_PATH}/${SRA_ID}_"

${TOOL_STAR}/STAR \
--runThreadN 8 \
--genomeDir ${GENOME_HOME}/${STAR_INDEX_PATH} \
--readFilesIn ${DATA}/${READ_1} ${DATA}/${READ_2} \
--outFileNamePrefix ${OUTPUT_FILE}