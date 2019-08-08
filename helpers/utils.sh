#!/usr/bin/env bash
#===============================================================================
#
#          FILE: utils.sh
#
#         USAGE: source utils.sh
#
#   DESCRIPTION: Useful functions which will be repetately used can be
#   imported from this file
#
#        AUTHOR: Rohit Suratekar
#  ORGANIZATION: IIMCB
#       CREATED: Thursday 08 August 2019 13:23
#      REVISION:  1
#===============================================================================

set -o nounset # Treat unset variables as an error

source ./constants/others.sh

log() {
  # Log the details of it is enabled
  if [[ ${LOG_ENABLED} -eq 1 ]]; then
    echo "RLOG $(date +%y-%m-%d" "%H:%M:%S) :$1" >>"${LOG_FILE}"
  fi
  # Display if display is enabled
  if [[ ${LOG_DISPLAY} -eq 1 ]]; then
    echo "RLOG $(date +%y-%m-%d" "%H:%M:%S) :$1"
  fi
}

record() {
  echo "RLOG $(date +%y-%m-%d" "%H:%M:%S) :$1" >>"${RECORD}"
}

check_file() {
  path_given="."
  while getopts ":p:" option; do
    case $option in
    p)
      echo "Path is given :  $OPTARG"
      ;;
    *)
      echo "None : $OPTARG"
      ;;
    esac
  done
  shift $((OPTIND - 1))

}
