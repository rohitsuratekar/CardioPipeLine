#!/usr/bin/env bash
#===============================================================================
#
#          FILE: salmon_post_processing.sh
#
#         USAGE: source ./salmon_post_processing.sh
#
#   DESCRIPTION: Post process the salmon output files for further analysis
#
#        AUTHOR: Rohit Suratekar
#  ORGANIZATION: IIMCB
#       CREATED: Wednesday 21 August 2019 11:28
#      REVISION:  1
#===============================================================================

set -o nounset # Treat unset variables as an error

source config.sh

python3 ./python/get_genes.py "${FOLDER_GENOME}/${GENOME_ANNOTATION}" \
  "${FOLDER_ANALYSIS}/${SRA_ID}/salmon/quant.sf"
