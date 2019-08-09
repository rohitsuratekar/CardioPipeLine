#!/usr/bin/env bash
#===============================================================================
#
#          FILE: star_indexing.sh
#
#         USAGE: ./star_indexing.sh
#
#   DESCRIPTION: Indexes Genome files for STAR alignment
#
#  REQUIREMENTS: Genome assembly and annotation files
#        AUTHOR: Rohit Suratekar
#  ORGANIZATION: IIMCB
#       CREATED: Friday 09 August 2019 17:16
#      REVISION:  1
#===============================================================================

set -o nounset # Treat unset variables as an error

source config.sh

perform_indexing() {
  log "Deleting all existing STAR indices"

  # Delete all existing files
  rm -r -f "${FOLDER_STAR_INDEX}"

  # Make new folder for keeping index files
  mkdir -p "${FOLDER_STAR_INDEX}"

  log "STAR Indexing started"

  "${TOOL_STAR}"/STAR \
    --runThreadN 8 \
    --runMode genomeGenerate \
    --genomeDir "${FOLDER_STAR_INDEX}" \
    --genomeFastaFiles "${FOLDER_GENOME}/${GENOME_FASTA}" \
    --sjdbGTFfile "${FOLDER_GENOME}/${GENOME_ANNOTATION}"
}

if perform_indexing ; then
  log "STAR Indexing completed"
else
  log "STAR indexing failed"
  exit 1
fi