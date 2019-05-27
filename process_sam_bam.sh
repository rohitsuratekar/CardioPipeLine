#!/usr/bin/env bash
#
# Processes sam files created by the alignment
# -T : prefix of the file
# -O : Type of output
# -o : output file. This will avoid printing file content to screen

source ./config.sh

readonly SAM_SUFFIX="Aligned.out.sam"
readonly INPUT_SAM_FILE="${DATA}/aligned/${SRA_ID}_${SAM_SUFFIX}"

log "Sorting of SAM file and conversion into BAM"

# Sort the sam files and convert to bam
${TOOL_SAMTOOL}/samtools sort \
-T "${DATA}/aligned/${SRA_ID}_sorted" \
-O bam \
-o "${DATA}/aligned/${SRA_ID}_sorted.bam" \
${INPUT_SAM_FILE}

log "Indexing of BAM file"
# Index the bam files

${TOOL_SAMTOOL}/samtools index "${DATA}/aligned/${SRA_ID}_sorted.bam"