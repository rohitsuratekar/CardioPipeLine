#!/usr/bin/env bash
#
# This script will download SRA file from the NCBI sever

source config.sh

# Download function for SRA tool
download_sra() {
  # Use prefetch from SRA tools to download file
  # -v : verbose output
  # -p 1: Show progress
  "${TOOL_SRA}"/bin/prefetch -v -p 1 "${SRA_ID}"
}

log "Preparing Files"
log "Checking if file is already downloaded"

# Check if file is in raw folder
if [[ -f "${DATA}/raw/${SRA_ID}.sra" ]]; then
  log "File '${SRA_ID}.sra' already exists in raw data. Skipping download."
  exit 0
# Check file is in NCBI's local downloading folder
elif [[ -f "${NCBI_DOWNLOAD}/${SRA_ID}.sra" ]]; then
  log "${SRA_ID}.sra file found in ncbs local folder. Moving to raw folder"
  mv "${NCBI_DOWNLOAD}/${SRA_ID}.sra" "${DATA}"/raw/
# If not, download with prefetch
else
  log "File '${SRA_ID}.sra' not found. Downloading with prefetch"
  download_sra
  log "${SRA_ID}.sra file downloaded"
  log "${SRA_ID}.sra file moving to raw folder"
  sleep 5
  mv "${NCBI_DOWNLOAD}/${SRA_ID}.sra" "${DATA}"/raw/
fi

# Check if script is successful
if [[ "$?" -eq 1 ]]; then
  log "Unable to download the file. Exiting script"
  exit 1
else
  log "File '${SRA_ID}.sra' is ready for further processing"
fi
