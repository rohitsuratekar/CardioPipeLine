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
#SRA_ID="SRX4720625"

# Check the quality of the sequence when it is raw and after it is filtered
CHECK_SEQ_QUALITY=1         # 0=NO, 1=YES

# Filters rRNA with SortMeRNA
RNA_FILTERING=0             # 0=NO, 1=YES

# Index salmon. If Yes, it will reindex the existing index (if any)
INDEX_SALMON=0              # 0=No, 1=YES

# Source All constants and extensions
source ./constants/paths.sh #Imports various paths
source ./constants/files.sh #Imports file names and extensions
source ./helpers/utils.sh   #Imports important functions
