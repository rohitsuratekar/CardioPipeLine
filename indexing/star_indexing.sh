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

# Delete all existing files
rm -r -f ${GENOME_HOME}/index/star

# Make new folder for keeping index files
mkdir -p ${GENOME_HOME}/index/star

readonly STAR_INDEX_PATH="index/star"
readonly STAR_REF_FILE="Danio_rerio.GRCz11.dna.only_chromosomes.fa"
readonly STAR_ANNOTATION_FILE="Danio_rerio.GRCz11.96.chr.gtf"


${TOOL_STAR}/STAR \
--runThreadN 8 \
--runMode genomeGenerate \
--genomeDir ${GENOME_HOME}/${STAR_INDEX_PATH} \
--genomeFastaFiles ${GENOME_HOME}/${STAR_REF_FILE} \
--sjdbGTFfile ${GENOME_HOME}/${STAR_ANNOTATION_FILE}