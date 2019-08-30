#!/bin/bash -
#===============================================================================
#
#          FILE: process_bam.sh
#
#         USAGE: ./process_bam.sh
#
#   DESCRIPTION: Uses StringTie to process the BAM files
#
#  REQUIREMENTS: Indexed and sorted BAM files
#         NOTES: You might need to adjust name of your BAM file if you have
#                processed it seperately from this pipeline
#        AUTHOR: Rohit Suratekar
#  ORGANIZATION: IIMCB
#       CREATED: Friday 09 August 2019 17:23
#      REVISION:  0
#===============================================================================

set -o nounset # Treat unset variables as an error

source config.sh

# Adjust this accordingly
BAM_FILENAME="${SRA_ID}_Aligned.sortedByCoord.out.bam"

stringtie_analysis() {

  # Get statistics with StringTie
  # -p : number of core
  # -G : Reference annotation
  # -e : only estimate the abundance of given reference transcripts (requires -G). You might want to ommit this if working with novel transcripts
  # -B : enable output of Ballgown table files
  # -o : output path/file name for the assembled transcripts GTF
  # -A : gene abundance estimates (provide output path/file name)

  log "StringTie analysis started for $SRA_ID"

  "${TOOL_STRINGTIE}/stringtie" \
    -p 8 \
    -G "${FOLDER_GENOME}/${GENOME_ANNOTATION}" \
    -e \
    -B \
    -o "${1}" \
    -A "${2}" \
    "${FOLDER_ALIGNED}/${SRA_ID}/${BAM_FILENAME}"

}

output_folder="${FOLDER_ANALYSIS}/${SRA_ID}"
mkdir -p "$output_folder"

GTF_OUTPUT_FILE="${output_folder}/${SRA_ID}_assembled.gtf"
GENE_EXPRESSION_OUTPUT="${output_folder}/${SRA_ID}_gene_expression.tsv"

if stringtie_analysis "$GTF_OUTPUT_FILE" "$GENE_EXPRESSION_OUTPUT"; then
  log "Analysis Complete, all output files are in ${output_folder}"
  log "Recommended: run clean up script for removing all intermediate files"
else
  log "StringTie analysis failed"
  exit 1
fi
