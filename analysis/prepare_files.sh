#!/usr/bin/env bash
#===============================================================================
#
#          FILE: prepare_files.sh
#
#         USAGE: ./prepare_files.sh
#
#   DESCRIPTION: Prepares downloaded SRA files for further processing
#
#       OPTIONS: ---
#  REQUIREMENTS: Downloaded SRA files in '$FOLDER_DATA/raw' folder
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Rohit Suratekar
#  ORGANIZATION: IIMCB
#       CREATED: Friday 09 August 2019 14:48
#      REVISION:  3
#===============================================================================

set -o nounset # Treat unset variables as an error

source config.sh

convert_to_fastq() {

  log "Starting SRA to FASTQ conversion"
  
  # Following command is borrowed from  https://edwards.sdsu.edu/research/fastq-dump/
  #
  #--gzip : Compression format (not required here)
  #--skip-technical : Exclude technical reads and includes only biological reads
  #--readids : Append IDs if fies are split by spot
  #--dumpbase : Formats sequence using base space (only required for SOLiD)
  #--split-files : Dump each read into separate file.
  #--clip : Removes left and right clips/tags
  #--outdir : If we want output in different dictionary
  #PATH/TO/SRA_FILE : Full path of file

  "${TOOL_SRA}/bin/fastq-dump" \
    --skip-technical \
    --readids \
    --dumpbase \
    --split-files \
    --clip \
    --outdir "${FOLDER_DATA}" \
    "${FOLDER_DATA}/raw/${SRA_ID}.sra"

}

# As there will be two file (forward and reverse) for our analysis, we will
# use file pattern to check if file exists
if check_file -p "$FOLDER_DATA" "${SRA_ID}*.fastq"; then
  log "FASTQ file(s) exists. Skipping FASTQ conversion."
  log "Check manually if you have 2 files in the data folder. If not,
  delete the existing files and run this script again"
  exit 0
fi

# If files are not there then use fastq-dump to convert them

if convert_to_fastq; then
  log "SRA files converted to FASTQ for further analysis"
  exit 0
else
  log "Unable to convert into FASTQ files"
  exit 1
fi
