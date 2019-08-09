#!/usr/bin/env bash
#===============================================================================
#
#          FILE: check_quality.sh
#
#         USAGE: ./check_quality.sh -r
#                ./check_quality.sh -f
#
#   DESCRIPTION: Checkes the quality of the FASTQ files provided
#
#       OPTIONS: -r (for raw) or -f (for filtered)
#  REQUIREMENTS: FASTQ files
#         NOTES: Can be skipped by setting $CHECK_SEQ_QUALITY to 0
#        AUTHOR: Rohit Suratekar
#  ORGANIZATION: IIMCB
#       CREATED: Friday 09 August 2019 15:05
#      REVISION:  3
#===============================================================================

set -o nounset # Treat unset variables as an error

source config.sh

# If it is configure to skip the sequence quality check, then just retuen
# the script as a successful

if [[ ${CHECK_SEQ_QUALITY} -eq 0 ]]; then
  log "Skipping sequence quality check. If you want to check sequence
  quality then change CHECK_SEQ_QUALITY value to 1 in config.sh file."
  exit 0
fi

# Check arguments

file_path="."
quality_suffix="."
unset OPTARG OPTIND
while getopts ":rf" option; do
  case $option in
  r)
    file_path="$FOLDER_DATA"
    quality_suffix="_RAW"
    ;;
  f)
    file_path="$FOLDER_FILTERING"
    quality_suffix="_FILTERED"
    ;;
  *)
    log "Unknown argument. use '-r' (for raw sequence) or '-f' (for
    filtered sequence)"
    exit 1
    ;;
  esac
done

# If no argument is passed, exit is script
if [[ "$file_path" == "." ]]; then
  log "For checking quality, passing argument is mandatory"
  log "use '-r' (for raw sequence) or '-f' (for filtered sequence)"
  exit 1
fi

check_seq_quality() {
  # Chekc quality with fastQC
  # --outdir : Folder for kepeeing the output files
  if "${TOOL_FASTQC}"/fastqc "${1}" --outdir="${2}"; then
    log "Quality checked for '${1}' ."
    log "results can be found in '${2}'"
  else
    log "Unable to check the quality of the sequence."
    exit 1
  fi
}

# Make folder to keep analysis files
output_folder="${FOLDER_QUALITY}/${SRA_ID}/${SRA_ID}${quality_suffix}"
mkdir -p "$output_folder"

# Check if fastq file exists
if ! check_file -p "${file_path}" "${SRA_ID}*.fastq"; then
  log "Unable to find the required FASTQ files. Please check if you have
  given correct argument"
  exit 1
fi

# Check quality
find "${file_path}" -maxdepth 1 -type f -name "${SRA_ID}*.fastq" |
  while read -r filename; do
    check_seq_quality "${filename}" "${output_folder}"
  done
