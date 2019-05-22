#!/usr/bin/env bash
# This script checks requirement of pipeline. 

source config.sh

# Check if all the tools are available

log "Checking availability of required tools"

${TOOL_SRA}/bin/fastq-dump --version
if ! [[ -e "${TOOL_HISAT}/hisat2_extract_splice_sites.py" ]]; then
     echo "hisat2_extract_splice_sites.py does not exist"
     exit 1
fi
if ! [[ -e "${TOOL_HISAT}/extract_exons.py" ]]; then
     echo "extract_exons.py does not exist"
     exit 1
fi
${TOOL_HISAT}/hisat2 --version
${TOOL_HISAT}/hisat2-build --version
${TOOL_FASTQC}/fastqc --version

if ! [[ -e "${TOOL_SORTMERNA}/scripts/merge-paired-reads.sh" ]]; then
     echo "scripts/merge-paired-reads.sh does not exist"
     exit 1
fi

if ! [[ -e "${TOOL_SORTMERNA}/scripts/unmerge-paired-reads.sh" ]]; then
     echo "scripts/merge-paired-reads.sh does not exist"
     exit 1
fi

if ! [[ -e "${TOOL_SORTMERNA}/bin/indexdb" ]]; then
     echo "bin/indexdb does not exist"
     exit 1
fi

if ! [[ -e "${TOOL_SORTMERNA}/bin/sortmerna" ]]; then
     echo "bin/indexdb does not exist"
     exit 1
fi

log "All Required Tools are installed"
