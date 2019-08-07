#!/usr/bin/env bash
#
# STAR indexing

source config.sh

readonly STAR_INDEX_PATH="index/star"

log "Deleting all existing STAR indices"

# Delete all existing files
rm -r -f "${GENOME_HOME:?}/${STAR_INDEX_PATH}"

# Make new folder for keeping index files
mkdir -p "${GENOME_HOME}/${STAR_INDEX_PATH}"

log "STAR Indexing started"

"${TOOL_STAR}"/STAR \
  --runThreadN 8 \
  --runMode genomeGenerate \
  --genomeDir "${GENOME_HOME}/${STAR_INDEX_PATH}" \
  --genomeFastaFiles "${GENOME_HOME}/${GENOME_FASTA}" \
  --sjdbGTFfile "${GENOME_HOME}/${GENOME_ANNOTATION}"

log "STAR Indexing completed"
