#!/usr/bin/env bash
#
# Use to delete all unnecessary file


source ./config.sh

log "Deleting unnecessary files from the analysis"

log "Deleting merged FASTQ file"
rm -rf "${DATA}/filtering/merged"

log "Deleting filtered merge files and rRNA files"
rm -rf "${DATA}/filtering/${SRA_ID}${R_RNA_FILTERED_EXTENSION}.fastq"
rm -rf "${DATA}/filtering/${SRA_ID}_rRNA.fastq"