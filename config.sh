#!/usr/bin/env bash
#===============================================================================
#
#          FILE: config.sh
#
#         USAGE: source config.sh
#
#   DESCRIPTION: All configuration settings can be set from this file
#
#        AUTHOR: Rohit Suratekar
#  ORGANIZATION: IIMCB
#       CREATED: Wednesday 07 August 2019 15:40
#      REVISION:  3
#===============================================================================

set -o nounset # Treat unset variables as an error

# SRA ID for the analysis, this ID will be used in whole Pipeline
SRA_ID="SRX4720634"

# Source All constants and extensions
source ./constants/paths.sh #Imports various paths
source ./constants/files.sh #Imports file names and extensions
source ./helpers/utils.sh   #Imports important functions
