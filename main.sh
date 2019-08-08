#!/usr/bin/env bash
#===============================================================================
#
#          FILE: main.sh
# 
#         USAGE: ./main.sh 
# 
#   DESCRIPTION: Main File to run the PipeLine
# 
#       OPTIONS: ---
#  REQUIREMENTS: ---
#          BUGS: ---
#         NOTES: ---
#        AUTHOR: Rohit Suratekar 
#  ORGANIZATION: IIMCB
#       CREATED: Wednesday 07 August 2019 15:20
#      REVISION: 2
#===============================================================================

set -o nounset                              # Treat unset variables as an error



bash ./analysis/download_files.sh
