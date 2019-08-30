#!/bin/bash -
#===============================================================================
#
#          FILE: clean_files.sh
#
#         USAGE: ./clean_files.sh
#
#   DESCRIPTION: Cleans Up all the files after analysis is over
#
#         NOTES: BE CAREFUL
#        AUTHOR: Rohit Suratekar
#  ORGANIZATION: IIMCB
#       CREATED: Friday 09 August 2019 17:37
#      REVISION:  1
#===============================================================================

set -o nounset # Treat unset variables as an error

source config.sh

log "Deleting unnecessary files from the analysis of ${SRA_ID}"

log "Deleting split FASTQ files"
rm -rf "${FOLDER_DATA}/${SRA_ID}_1.fastq"
rm -rf "${FOLDER_DATA}/${SRA_ID}_2.fastq"

log "Deleting aligned files"
rm -rf "${FOLDER_DATA}/aligned/${SRA_ID}"

log "Deleting filtering data"
rm -rf "${FOLDER_FILTERING}"

log "Deleting downloaded files"
rm -rf "${FOLDER_DATA}/raw/${SRA_ID}.sra"
