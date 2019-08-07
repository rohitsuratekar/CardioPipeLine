#!/usr/bin/env bash
#
# Author : Rohit Suratekar
# Date : July 2019
#
# Functions related to FASTQC quality check
# Suffix arguments can be given while running the script
#
# -r = raw files
# -f = filtered files
#
# e.g : check_quality.sh -f

source ./config.sh

# If it is configure to skip the sequence quality check, then just retuen
# the script as a successful

if [[ ${CHECK_SEQ_QUALITY} -eq 0 ]]; then
  log "Skipping sequence quality check. If you want to check sequence
  quality then change CHECK_SEQ_QUALITY value to 1 in config.sh file."
  exit 0
fi

check_seq_quality() {
  # Chekc quality with fastQC
  # --outdir : Folder for kepeeing the output files

  "${TOOL_FASTQC}"/fastqc "${1}" --outdir="${2}"
}

QUALITY_SUFFIX="_RAW"
QUALITY_INPUT="${DATA}"
# Check for any argument given and make the prefix out of it
if [[ $1 == "-f" ]]; then
  QUALITY_SUFFIX="_FILTERED"
  QUALITY_INPUT="${FILTERED_FOLDER}"
fi

# Make folder to keep analysis files
QUALITY_FOLDER="${DATA}/quality/${SRA_ID}/${SRA_ID}${QUALITY_SUFFIX}"
mkdir -p "${QUALITY_FOLDER}"

find "${QUALITY_INPUT}" -maxdepth 1 -name "${SRA_ID}*.fastq" -size +1M >/dev/null
if [[ "$?" -eq 1 ]]; then
  log "Unable to calculate the quality. Please check if '-f' argument is
  given in case of filtered files. Exiting Script."
  exit 1
fi

# Check quality
find "${QUALITY_INPUT}" -maxdepth 1 -name "${SRA_ID}*.fastq" -size +1M |
  while read -r filename; do
    check_seq_quality "${filename}" "${QUALITY_FOLDER}"
  done

# Check if script is successful
if [[ "$?" -ne 0 ]]; then
  log "Some problem in checking sequence quality. Exiting script"
  exit 1
else
  log "Quality checked for FASTQ files. results can be found in '${QUALITY_FOLDER}'"
fi
