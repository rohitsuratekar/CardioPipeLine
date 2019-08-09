#!/usr/bin/env bash
#===============================================================================
#
#          FILE: others.sh
#
#         USAGE: others.sh
#
#   DESCRIPTION: All other constants will go here
#
#        AUTHOR: Rohit Suratekar
#  ORGANIZATION: IIMCB
#       CREATED: Thursday 08 August 2019 12:20
#      REVISION: 3
#===============================================================================

set -o nounset            # Treat unset variables as an error

# Log related options
LOG_FILE="script_log.txt" # Name of the log file
LOG_ENABLED=0             # 1=Yes and 0=NO (Add to Log File)
LOG_DISPLAY=1             # 1=Yes and 0=No (Display on Terminal)

SLEEP_TIME=1 # Delay between two scripts. Used in few scripts
