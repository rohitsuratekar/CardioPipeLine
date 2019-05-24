#!/usr/bin/env bash
#
# Make Alignment of RNA-seq paired end data to the reference genome index using HISAT2
#
# Command was adapted from :
# https://github.com/griffithlab/rnaseq_tutorial/wiki/Alignment
# and
# https://ccb.jhu.edu/software/hisat2/manual.shtml
#
# -p : Number of CPU Cores
# -f :
# --rg-id : Unique Identifier for reads (will be added to BAM file)
# --rg SM : Sample name (will be added to BAM file)
# --rg PL : Sequencing Platform (will be added to BAM file)
# -x HISAT index
# -dta : Reports alignments tailored for transcript assemblers.
# --rna-strandness Rf : specifies strandness of RNAseq library. We will specify RF since the TruSeq strand-specific library was used to make these libraries
# -1 : Read 1
# -2 : Read 2
# -S : Output SAM file

source ./config.sh

readonly INDEX_FILE="zebrafish_index.fasta"
readonly SAMPLE_NAME="control_GFP"
readonly READ_1="SRX4720625_pass_1.fastq"
readonly READ_2="SRX4720625_pass_2.fastq"
readonly OUTPUT_ALIGN="GFP_POS1"
readonly SEQ_PLATFORM="ILLUMINA"

#${TOOL_HISAT}/hisat2 \
#-p 8 \
#--rg-id=${OUTPUT_ALIGN} \
#--rg SM:${SAMPLE_NAME} \
#--rg PL:${SEQ_PLATFORM} \
#-x ${GENOME_HOME}/index/${INDEX_FILE} \
#--dta \
#--rna-strandness RF \
#-1 ${DATA}/${READ_1} \
#-2 ${DATA}/${READ_2} \
#-S ${DATA}/${OUTPUT_ALIGN}.sam


${TOOL_HISAT}/hisat2 \
-p 8 \
-x ${GENOME_HOME}/index/${INDEX_FILE} \
-1 ${DATA}/${READ_1} \
-2 ${DATA}/${READ_2} \
-S ${DATA}/${OUTPUT_ALIGN}.sam