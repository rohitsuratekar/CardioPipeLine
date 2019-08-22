#!/usr/bin/env bash
#===============================================================================
#
#          FILE: salmon_post_processing.sh
#
#         USAGE: source ./salmon_post_processing.sh
#
#   DESCRIPTION: Post process the salmon output files for further analysis
#
#        AUTHOR: Rohit Suratekar
#  ORGANIZATION: IIMCB
#       CREATED: Wednesday 21 August 2019 11:28
#      REVISION:  1
#===============================================================================

set -o nounset # Treat unset variables as an error

source config.sh

INPUT_FILE="${FOLDER_ANALYSIS}/${SRA_ID}/salmon/quant.sf"
OUTPUT_FILE="${FOLDER_ANALYSIS}/${SRA_ID}/salmon/converted_${SRA_ID}.tsv"

run_python() {
  python3 ./external/get_genes.py "${FOLDER_GENOME}/${TRANSCRIPTS_TO_GENE}" \
    "$INPUT_FILE" "$OUTPUT_FILE"
}

log "Starting post-processing on the salmon file"

if run_python ; then
  log "Processed salmon file $OUTPUT_FILE"
else
  log "Could not process the salmon output file"
  return 1
fi