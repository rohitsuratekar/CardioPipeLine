#!/usr/bin/env bash
#===============================================================================
#
#          FILE: download_files.sh
# 
#         USAGE: ./download_files.sh 
# 
#   DESCRIPTION: 
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Rohit Suratekar 
#  ORGANIZATION: IIMCB
#       CREATED: Thursday 08 August 2019 12:20
#      REVISION: 3
#===============================================================================

set -o nounset                              # Treat unset variables as an error


source config.sh

# Download function for SRA tool
download_sra() {
  # Use prefetch from SRA tools to download file
  # -v : verbose output
  # -p 1: Show progress
  "${TOOL_SRA}"/bin/prefetch -v -p 1 "${SRA_ID}"
}


check_file "Name" -p "AOmr:path"

