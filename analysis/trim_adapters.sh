#!/usr/bin/env bash
#
# Trim adapters and improve the quality of the sequence

source config.sh

# PE: PairEnd
# -threads : Number of threads

FWD_SEQ="${DATA}/${SRA_ID}_1.fastq"
REV_SEQ="${DATA}/${SRA_ID}_2.fastq"
ADAPTER_FILE="${TOOL_TRIMMOMATIC}/adapters/TruSeq3-PE.fa"
TRIMMED_PREFIX="${DATA}/${SRA_ID}_TRIMED_"

java -jar "${TOOL_TRIMMOMATIC}/trimmomatic-0.39.jar" \
  PE \
  -threads 8 \
  "${FWD_SEQ}" "${REV_SEQ}" \
  -baseout "${TRIMMED_PREFIX}" \
  ILLUMINACLIP:"${ADAPTER_FILE}":2:30:10 \
  LEADING:3 \
  TRAILING:3 \
  SLIDINGWINDOW:4:15 \
  MINLEN:36

# Check if script is successful
if [[ "$?" -ne 0 ]]; then
  log "Unable to trim sequence files. Exiting script"
  exit 1
else
  log "Sequence files trimmed for further analysis"
fi
