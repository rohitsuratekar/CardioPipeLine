#!/usr/bin/env bash
#===============================================================================
#
#          FILE: main.sh
#
#         USAGE: ./main.sh
#
#   DESCRIPTION: Main File to run the PipeLine
#
#         NOTES: General pipline for RNA-seq analysis
#        AUTHOR: Rohit Suratekar
#  ORGANIZATION: IIMCB
#       CREATED: Wednesday 07 August 2019 15:20
#      REVISION: 3
#===============================================================================

stop() {
  echo "Something went wrong..."
  echo "Exiting script"
  exit 1
}
trap 'stop' ERR

arrange_files() {
  bash ./analysis/download_files.sh
  bash ./analysis/prepare_files.sh
  bash ./analysis/check_quality.sh -r
  bash ./analysis/index_rrna.sh
  bash ./analysis/rrna_filtering.sh
  bash ./analysis/check_quality.sh -f
}

star() {
  bash ./analysis/star_alignment.sh
  bash ./analysis/process_bam.sh
}

clean() {
  bash ./analysis/clean_files.sh
}

salmon() {
  bash ./analysis/salmon.sh
  bash ./analysis/salmon_post_processing.sh
}

pipeline() {
  while IFS='' read -r line; do
    current_id="$(echo -e "${line}" | tr -d '[:space:]')"
    export SRA_ID="$current_id"
    arrange_files
    star
    salmon
    clean
  done <"sra_ids.txt"
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
  pipeline
fi
