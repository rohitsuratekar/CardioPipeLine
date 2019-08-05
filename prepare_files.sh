#!/usr/bin/env bash
#
# Converts SRA files to FASTQ for further analysis

source config.sh

convert_to_fastq() {

  # Following command is borrowed from  https://edwards.sdsu.edu/research/fastq-dump/
  #
  #--gzip : Compression format (not required here)
  #--skip-technical : Exclude technical reads and includes only biological reads
  #--readids : Append IDs if fies are split by spot
  #--dumpbase : Formats sequence using base space (only required for SOLiD)
  #--split-files : Dump each read into separate file.
  #--clip : Removes left and right clips/tags
  #--outdir : If we want output in different dictionary
  #PATH/TO/SRA_FILE : Full path of file

  "${TOOL_SRA}/bin/fastq-dump" \
    --skip-technical \
    --readids \
    --dumpbase \
    --split-files \
    --clip \
    --outdir "${DATA}" \
    "${DATA}/raw/${SRA_ID}.sra"

}

log "Checking if FASTQ file exists"

# -maxdepth : look only in current folder level
# -name : exact name
# -size : only files greater than 1MB
# wc -l : count lines coming from find output

readonly TEMP1=$(find "${DATA}" \
  -maxdepth 1 \
  -name "${SRA_ID}*.fastq" \
  -size +1M |
  wc -l)

if [[ ${TEMP1} -gt 0 ]]; then
  log "Total of ${TEMP1} FASTQ files exists. Skipping FASTQ conversion."
  exit 0
else
  log "FASTQ files not found. Starting SRA to FASTQ conversion"
  convert_to_fastq
fi

# Check if script is successful
if [[ "$?" -eq 1 ]]; then
  log "Unable to convert into FASTQ files. Exiting script"
  exit 1
else
  log "SRA files converted to FASTQ for further analysis"
fi
