#!/usr/bin/env bash
#
# Author : Rohit Suratekar
# Date : July 2019
#
# All local level configuration can be adjusted from here. File names and
# extensions can be adjusted from files from 'constants' folder

# Following will import constant called LOCAL_FOLDER which is essentially
# path to your working dictionary. something like /home/user/Documents/work
source ./helpers/logging.sh
source ./constants/paths.sh


# SRA ID for the analysis
SRA_ID="SRX4720634"








# All log related constants can be edited from the helpers/logging.sh file
#
#log "===================================="
#log "Script: $(basename "$0")"
#log "Time: $(date)"
#log "SRA ID: ${SRA_ID}"
#log "===================================="
