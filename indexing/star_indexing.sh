#!/usr/bin/env bash
#
# Generates STAR index
#
# --runThreadN : Number of threads
# --runMode : genomeGenerate (Mode to generate indices)
# --genomeDir : Name of output dictionary
# --genomeFastaFiles : Reference Genome Fasta files
# --sjdbGTFfile : Annotation file
#
# Note: IDE might give warning for config.sh file but you can ignore them.

source ./config.sh

readonly STAR_INDEX_PATH="index/star"

log "Deleting all existing STAR indices"

# Delete all existing files
rm -r -f "${GENOME_HOME}/${STAR_INDEX_PATH}"

# Make new folder for keeping index files
mkdir -p "${GENOME_HOME}/${STAR_INDEX_PATH}"


log "STAR Indexing started"


${TOOL_STAR}/STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir ${GENOME_HOME}/${STAR_INDEX_PATH} \
--genomeFastaFiles ${GENOME_HOME}/${GENOME_FASTA} \
--sjdbGTFfile ${GENOME_HOME}/${GENOME_ANNOTATION}


log "STAR Indexing completed"