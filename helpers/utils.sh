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
  #Log function
  # To set these values use constants/other.sh file

  # Log the details of it is enabled
  if [[ ${LOG_ENABLED} -eq 1 ]]; then
    echo "RLOG $(date +%y-%m-%d" "%H:%M:%S) :$1" >>"${LOG_FILE}"
  fi
  # Display if display is enabled
  if [[ ${LOG_DISPLAY} -eq 1 ]]; then
    echo "RLOG $(date +%y-%m-%d" "%H:%M:%S) :$1"
  fi
}

check_file() {
  # Checks if file exists in given path and returns the status code
  #
  # -p = Path (default : . )
  # -s = Size of the file (default : +1M )

  # Reset the OPTARG and OPTIND to use getopts
  unset OPTARG OPTIND

  path_given="."
  size="+1M"
  while getopts ":p:s:" option; do
    case $option in
    p)
      path_given="$OPTARG"
      ;;
    s)
      size="$OPTARG"
      ;;
    *)
      echo "Unknown argument. use '-p' for the path"
      ;;
    esac
  done
  shift $((OPTIND - 1))

  file_name="$*"

  # -maxdepth : look only in current folder level
  # -name : exact name or pattern
  # -size : only files greater than 1MB
  # -type : Check only files
  # wc -l : count lines coming from find output

  n=$(find "$path_given" -maxdepth 1 -type f -name "$file_name" -size "$size" | wc -l)

  if [[ "$n" -gt 0 ]]; then
    return 0
  else
    return 1
  fi

}
