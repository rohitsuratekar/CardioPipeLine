#!/usr/bin/env bash
#
# Author : Rohit Suratekar
# Date : July 2019
#
# All local level configuration can be adjusted from here. File names and
# extensions can be adjusted from files from 'constants' folder

source ./helpers/logging.sh #Imports Logging functions
source ./constants/paths.sh #Imports various paths
source ./constants/files.sh #Imports file names and extensions

# SRA ID for the analysis
SRA_ID="SRX4720634"

CHECK_SEQ_QUALITY=0 # 0=NO, 1 = YESv
RNA_FILTERING=1     # 1=Yes and 0=No (If Yes, it will perform rRNA Filtering)

SCRIPT_DELAY=2 # Delay between two scripts (in seconds)

# All log related constants can be edited from the helpers/logging.sh file
#
#log "===================================="
#log "Script: $(basename "$0")"
#log "Time: $(date)"
#log "SRA ID: ${SRA_ID}"
#log "===================================="
