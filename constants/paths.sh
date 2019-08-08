#!/usr/bin/env bash
#===============================================================================
#
#          FILE: paths.sh
#
#         USAGE: source paths.sh
#
#   DESCRIPTION: All path constants can be adjusted from this file
#
#        AUTHOR: Rohit Suratekar
#  ORGANIZATION: IIMCB
#       CREATED: Wednesday 07 August 2019 15:40
#      REVISION:  3
#===============================================================================

source local.sh

# Various system paths
ORIGIN=${LOCAL_FOLDER}                      # Base folder for everything e.g /home/user/path/to/base/folder
FOLDER_GENOME="${ORIGIN}/Genome"            # Standard genomes and their indexes
TOOL="${ORIGIN}/Tools"                      # Tools involved in this pipeline
FOLDER_PIPE="${ORIGIN}/Pipeline"            # All pipeline scripts
FOLDER_DATA="${ORIGIN}/Data"                # All experimental data files
FOLDER_NCBI_DOWNLOAD=${NCBI_LOCAL}          # /home/user/ncbi/public/sra

# Various analysis related paths
FOLDER_RAW_DATA="${FOLDER_DATA}/raw"        # Here all downloaded raw data will be stored
FOLDER_FILTERING="${FOLDER_DATA}/filtering" # To keep all filtered files
FOLDER_ANALYSIS="${FOLDER_DATA}/analysis"   # TO Keep final analysis files

# All Tools paths
TOOL_SRA="$TOOL/sratoolkit.2.9.6-1-ubuntu64"
TOOL_HISAT="$TOOL/hisat2-2.1.0"
TOOL_HT_SEQ="$TOOL/htseq-release_0.11.0"
TOOL_SORTMERNA="$TOOL/sortmerna-3.0.3"
TOOL_FASTQC="$TOOL/FastQC"
TOOL_STAR="$TOOL/STAR-2.7.1a/source"
TOOL_SAMTOOL="$TOOL/samtools-1.9"
TOOL_STRINGTIE="$TOOL/stringtie-1.3.4d.Linux_x86_64"
TOOL_TRIMMOMATIC="$TOOL/Trimmomatic-0.39"
