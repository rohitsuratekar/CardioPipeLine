#!/usr/bin/env bash
#
# Alignment and Mapping using STAR
#

source ./config.sh

star_mapping() {

  # --runThreadN : Number of threads to run
  # --genomeDir : Path of folder where indices are generated. (In this pipeline from "star_indexing.sh" script)
  # --readFilesIn : Reads to be mapped
  # --outFileNamePrefix : Output file path
  # --outSAMtype BAM SortedByCoordinate : Directly makes sorted BAM file so
  # there is no need to use samtools after this

  "${TOOL_STAR}"/STAR \
    --runThreadN 8 \
    --genomeDir "${GENOME_HOME}/${1}" \
    --readFilesIn "${2}" "${3}" \
    --outFileNamePrefix "${4}" \
    --outSAMtype BAM SortedByCoordinate

  # Check if command is successful
  if [[ "$?" -ne 0 ]]; then
    log "Unable to process STAR alignment. Exiting script"
    exit 1
  else
    log "STAR Alignment finished"
  fi

}

STAR_INDEX_PATH="index/star"
READ_1="${FILTERED_FOLDER}/${SRA_ID}${R_RNA_FILTERED_EXTENSION}_1.fastq"
READ_2="${FILTERED_FOLDER}/${SRA_ID}${R_RNA_FILTERED_EXTENSION}_2.fastq"
OUTPUT_PATH="${DATA}/aligned/${SRA_ID}"

mkdir -p "${OUTPUT_PATH}"
OUTPUT_FILE="${OUTPUT_PATH}/${SRA_ID}_"

if [[ ${RNA_FILTERING} -eq 0 ]]; then
  log "rRNA filtering is not done. Using file without rRNA filtering"
  READ_1="${DATA}/${SRA_ID}_1.fastq"
  READ_2="${DATA}/${SRA_ID}_2.fastq"
else
  log "Initializing STAR Alignment"
fi

# Start the alignment
star_mapping "${STAR_INDEX_PATH}" "${READ_1}" "${READ_2}" "${OUTPUT_FILE}"

ANALYSIS_OUPUT="${FOLDER_FINAL_OUTPUT}/${SRA_ID}"
STAR_OUTPUT_LOG="${OUTPUT_PATH}/${SRA_ID}_Log.final.out"
STAR_SJ_FILE="${OUTPUT_PATH}/${SRA_ID}_SJ.out.tab"

# Copy files to the analysis folder
mkdir -p "${ANALYSIS_OUPUT}"
cp "${STAR_OUTPUT_LOG}" "${ANALYSIS_OUPUT}/"
cp "${STAR_SJ_FILE}" "${ANALYSIS_OUPUT}/"

# Check if script is successful
if [[ "$?" -ne 0 ]]; then
  log "Unable to copy STAR output files. Exiting script"
  exit 1
else
  log "STAR analysis files copied to ${ANALYSIS_OUPUT}"
fi
