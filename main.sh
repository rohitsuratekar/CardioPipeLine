#!/usr/bin/env bash
#===============================================================================
#
#          FILE: main.sh
#
#         USAGE: ./main.sh
#
#   DESCRIPTION: Main File to run the PipeLine
#
#         NOTES: ---
#        AUTHOR: Rohit Suratekar
#  ORGANIZATION: IIMCB
#       CREATED: Wednesday 07 August 2019 15:20
#      REVISION: 2
#===============================================================================

stop() {
  echo "Something went wrong..."
  echo "Exiting script"
  exit 1
}
trap 'stop' ERR

star() {
  bash ./analysis/download_files.sh
  bash ./analysis/prepare_files.sh
  bash ./analysis/check_quality.sh -r
  bash ./analysis/index_rrna.sh
  bash ./analysis/rrna_filtering.sh
  bash ./analysis/check_quality.sh -f
  bash ./analysis/star_alignment.sh
  bash ./analysis/process_bam.sh
}

clean() {
  bash ./analysis/clean_files.sh
}

salmon() {
  bash ./analysis/download_files.sh
  bash ./analysis/prepare_files.sh
  bash ./analysis/salmon.sh
}

arg=0
while getopts ":c" option; do
  case $option in
  c)
    arg=1
    ;;
  *)
    log "Unknown argument. use '-c' (for cleaning files)"
    exit 1
    ;;
  esac
done

if [[ "$arg" -eq 1 ]]; then
  clean
else
  salmon
fi
