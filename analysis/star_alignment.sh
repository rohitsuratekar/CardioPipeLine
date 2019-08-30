#!/usr/bin/env bash
#===============================================================================
#
#          FILE: star_alignment.sh
#
#         USAGE: ./star_alignment.sh
#
#   DESCRIPTION: Performs STAR alignment
#
#  REQUIREMENTS: STAR indexing
#        AUTHOR: Rohit Suratekar
#  ORGANIZATION: IIMCB
#       CREATED: Friday 09 August 2019 16:59
#      REVISION:  3
#===============================================================================

set -o nounset # Treat unset variables as an error

source config.sh

star_mapping() {

  # --runThreadN : Number of threads to run
  # --genomeDir : Path of folder where indices are generated. (In this pipeline from "star_indexing.sh" script)
  # --readFilesIn : Reads to be mapped
  # --outFileNamePrefix : Output file path
  # --outSAMtype BAM SortedByCoordinate : Directly makes sorted BAM file so
  # there is no need to use samtools after this

  log "START alignment started for $SRA_ID"

  "${TOOL_STAR}"/STAR \
    --runThreadN 8 \
    --genomeDir "$FOLDER_STAR_INDEX" \
    --readFilesIn "${1}" "${2}" \
    --outFileNamePrefix "${3}" \
    --outSAMtype BAM SortedByCoordinate

}

READ_1="${FOLDER_FILTERING}/${SRA_ID}${EXTENSION_R_RNA_FILTERED}_1.fastq"
READ_2="${FOLDER_FILTERING}/${SRA_ID}${EXTENSION_R_RNA_FILTERED}_2.fastq"

# if rRNA filtering if OFF, use unfilterd sequence for alignment
if [[ ${RNA_FILTERING} -eq 0 ]]; then
  log "rRNA filtering is not done."
  log "Using unfilterd sequence for STAR alignment"
  READ_1="${FOLDER_DATA}/${SRA_ID}_1.fastq"
  READ_2="${FOLDER_DATA}/${SRA_ID}_2.fastq"
fi

output_path="${FOLDER_ALIGNED}/${SRA_ID}"
mkdir -p "${output_path}"

output_file="${output_path}/${SRA_ID}_"

if star_mapping "${READ_1}" "${READ_2}" "${output_file}"; then
  log "STAR Alignment finished"
else
  log "Unable to process STAR alignment."
  exit 1
fi

# Copy important log files to the analysis folder

ANALYSIS_OUPUT="${FOLDER_ANALYSIS}/${SRA_ID}/mapping"
STAR_OUTPUT_LOG="${output_path}/${SRA_ID}_Log.final.out"
STAR_SJ_FILE="${output_path}/${SRA_ID}_SJ.out.tab"

# Copy files to the analysis folder
mkdir -p "${ANALYSIS_OUPUT}"
cp "${STAR_OUTPUT_LOG}" "${ANALYSIS_OUPUT}/"
cp "${STAR_SJ_FILE}" "${ANALYSIS_OUPUT}/"

log "STAR Alignment log files copied to $ANALYSIS_OUPUT"
