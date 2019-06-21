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

log "Deleting split FASTQ files"
rm -rf "${DATA}/${SRA_ID}_1.fastq"
rm -rf "${DATA}/${SRA_ID}_2.fastq"

log "Deleting aligned files"
rm -rf "${DATA}/aligned/${SRA_ID}"

log "Deleting filtering data"
rm -rf "${DATA}/filtering"