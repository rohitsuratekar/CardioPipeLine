#!/usr/bin/env bash
#===============================================================================
#
#          FILE: download_files.sh
#
#         USAGE: ./download_files.sh
#
#   DESCRIPTION: Downloads given SRA from the NCBI server
#
#        AUTHOR: Rohit Suratekar
#  ORGANIZATION: IIMCB
#       CREATED: Thursday 08 August 2019 12:20
#      REVISION: 3
#===============================================================================

set -o nounset # Treat unset variables as an error

source config.sh

# Download function for SRA tool
download_sra() {
  # Use prefetch from SRA tools to download file
  # -v : verbose output
  # -p 1: Show progress

  if "${TOOL_SRA}/bin/prefetch" -v -p 1 "$SRA_ID"; then
    log "${SRA_ID}.sra file downloaded"
    return 0
  else
    log "Unable to fetch the file from NCBI server."
    exit 1
  fi
}

is_downloaded() {
  # Checks if given SRA file is already downloaded
  # If file is available in the NCBI public folder, move it to the raw folder

  log "Checking if file is already downloaded"
  if check_file -p "$FOLDER_RAW_DATA" "${SRA_ID}.sra"; then
    log "File '${SRA_ID}.sra' already exists in raw data. Skipping download."
    exit 0
  elif check_file -p "$FOLDER_NCBI_DOWNLOAD" "${SRA_ID}.sra"; then
    log "${SRA_ID}.sra file found in ncbs local folder. Moving to raw folder"
    mv "${FOLDER_NCBI_DOWNLOAD}/${SRA_ID}.sra" "$FOLDER_RAW_DATA"
    log "${SRA_ID}.sra file moved to raw folder"
    exit 0
  else
    return 1
  fi
}

# First check if file already exists
# If it is present then skip download
if is_downloaded; then
  exit 0
else
  log "File '${SRA_ID}.sra' not found. Downloading with prefetch"
  if download_sra; then
    sleep "$SLEEP_TIME"
    log "${SRA_ID}.sra file moving to raw folder"
    mv "${FOLDER_NCBI_DOWNLOAD}/${SRA_ID}.sra" "$FOLDER_RAW_DATA"
    log "${SRA_ID}.sra file moved to raw folder"
    exit 0
  fi
fi
