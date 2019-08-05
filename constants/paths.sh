#!/usr/bin/env bash
#
# All the paths in current folders

source ./local.sh

# Various paths
ORIGIN=${LOCAL_FOLDER} # Base folder for everything e.g /home/user/path/to/base/folder
GENOME_HOME="$ORIGIN/Genome" # Standard genomes and their indexes
TOOL="$ORIGIN/Tools" # Tools involved in this pipeline
PIPE="$ORIGIN/Pipeline" # All pipeline scripts
DATA="$ORIGIN/Data" # All experimental data files
NCBI_DOWNLOAD=${NCBI_LOCAL} # /home/user/ncbi/public/sra
RAW_DATA="$DATA/raw" # Here all downloaded raw data will be stored

# All Tools paths

TOOL_SRA="$TOOL/sratoolkit.2.9.6-1-ubuntu64"
TOOL_HISAT="$TOOL/hisat2-2.1.0"
TOOL_HT_SEQ="$TOOL/htseq-release_0.11.0"
TOOL_SORTMERNA="$TOOL/sortmerna-3.0.3"
TOOL_FASTQC="$TOOL/FastQC"
TOOL_STAR="$TOOL/STAR-2.7.1a/source"
TOOL_SAMTOOL="$TOOL/samtools-1.9"
TOOL_STRINGTIE="$TOOL/stringtie-1.3.4d.Linux_x86_64"

# Additional requirements


mkdir -p "${RAW_DATA}" # Make folder for raw data if not present