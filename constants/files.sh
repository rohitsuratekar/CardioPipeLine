#!/usr/bin/env bash
#===============================================================================
#
#          FILE: files.sh
#
#         USAGE: source files.sh
#
#   DESCRIPTION: All constants related to file names and extension will go
#   in this files.
#
#        AUTHOR: Rohit Suratekar
#  ORGANIZATION: IIMCB
#       CREATED: Thursday 08 August 2019 10:00
#      REVISION:  3
#===============================================================================


# Genome related files

GENOME_FASTA="Danio_rerio.GRCz11.dna.primary_assembly.fa" #Name of the fasta reference file
GENOME_ANNOTATION="Danio_rerio.GRCz11.96.chr.gtf"         # Genome Annotation file

# File extensions

 # Extension which will be added at the end of rRNA index file
EXTENSION_R_RNA_INDEX=".idx"

# Extension which will be added at the end of aligned rRNA filtered sequence
EXTENSION_R_RNA_FILTERED=".filtered"

# Extension will be added at the end of merged file
EXTENSION_MERGED_READS="_merged.fastq"
