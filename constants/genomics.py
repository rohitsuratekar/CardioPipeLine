#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# All genomics files
# Please update files according to your need
#
# Two constants DEX_HOME, DEX_WORK and DEX_DATA are taken from local.py files
# which are private. You should change these paths here according to your
# system.
#
# Avoid trailing slash in the path. It should be handled at the usage level

from local import *

# This file has all the accessions for SRA database and their relations
# URL: https://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/
# Accession date: 09 December 2019
SRA_ACCESSION = f"{DEX_DATA}/annotations/SRA_Accessions.tab"

ZEBRAFISH_GENOME_FASTA = f"{DEX_WORK}/Genome/Danio_rerio.GRCz11.dna.primary_assembly.fa"
ZEBRAFISH_GTF_ANNOTATION = f"{DEX_WORK}/Genome/Danio_rerio.GRCz11.98.chr.gtf"
RRNA_DATABASE = f"{DEX_WORK}/Genome/rrna"  # Need for rrna filtering
