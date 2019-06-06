#!/usr/bin/env bash
#
# Processes sam files created by the alignment
# -T : prefix of the file
# -O : Type of output
# -o : output file. This will avoid printing file content to screen

source ./config.sh

readonly SAM_SUFFIX="Aligned.out.sam"
readonly INPUT_SAM_FILE="${DATA}/aligned/${SRA_ID}/${SRA_ID}_${SAM_SUFFIX}"

log "Sorting of SAM file and conversion into BAM"

# Sort the sam files and convert to bam
${TOOL_SAMTOOL}/samtools sort \
-T "${DATA}/aligned/${SRA_ID}/${SRA_ID}_sorted" \
-O bam \
-o "${DATA}/aligned/${SRA_ID}/${SRA_ID}_sorted.bam" \
${INPUT_SAM_FILE}

# Index the bam files
log "Indexing of BAM file"
${TOOL_SAMTOOL}/samtools index "${DATA}/aligned/${SRA_ID}/${SRA_ID}_sorted.bam"


# Generate FASTA reference index
log "Generating fasta index"

${TOOL_SAMTOOL}/samtools faidx "${GENOME_HOME}/${GENOME_FASTA}"


# Get statistics with StringTie
# -p : number of core
# -G : Reference annotation
# -e : only estimate the abundance of given reference transcripts (requires -G). You might want to ommit this if working with novel transcripts
# -B : enable output of Ballgown table files
# -o : output path/file name for the assembled transcripts GTF
# -A : output path/file name for gene abundance estimates

log "Generating Statistics"

readonly GTF_OUTPUT_PATH="${DATA}/analysis/${SRA_ID}"
mkdir -p ${GTF_OUTPUT_PATH}
readonly GTF_OUTPUT_FILE="${GTF_OUTPUT_PATH}/${SRA_ID}_assembled.gtf"
readonly GENE_EXPRESSION_OUTPUT="${GTF_OUTPUT_PATH}/${SRA_ID}_gene_expression.tsv"

${TOOL_STRINGTIE}/stringtie \
-p 8 \
-G ${GENOME_HOME}/${GENOME_ANNOTATION} \
-e \
-B \
-o ${GTF_OUTPUT_FILE} \
-A ${GENE_EXPRESSION_OUTPUT} \
${DATA}/aligned/${SRA_ID}/${SRA_ID}_sorted.bam


log "Analysis Complete, all output files are in ${GTF_OUTPUT_PATH}"
log "Recommended: run clean up script for removing all intermediate files"