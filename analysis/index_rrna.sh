#!/usr/bin/env bash
#
# Index rRNA required for the rRNA filtering

source config.sh

# If it is configure to skip the rRNA filtering, then just retuen
# the script as a successful

if [[ ${RNA_FILTERING} -eq 0 ]]; then
  log "Skipping rRNA indexing. If you want to filter sequences with rRNA
  then change RNA_FILTERING value to 1 in config.sh file."
  exit 0
fi

# Path of folder where you want your fasta files and index files to be stored
readonly R_RNA_PATH="$GENOME_HOME/rrna"

# Index files will be stored here
mkdir -p "${R_RNA_PATH}"/index

index_rrna() {
  # Following command is VERY sensitive to the spaces. So double check if
  # you are not giving any additional spaces in the command

  "${TOOL_SORTMERNA}/bin/indexdb" --ref "${1}","${R_RNA_PATH}/index/${2}${R_RNA_INDEX_EXTENSION}" -v # Make Index
  # Check if indexing is successful
  if [[ "$?" -ne 0 ]]; then
    log "Some problem in indexing rRNA. Exiting script"
    exit 1
  fi

}

# Check if rRNA FASTA files exists

readonly TEMP2=$(find "${R_RNA_PATH}" -name "*.fasta" | wc -l)

# This will skip copying files from SortMeRNA even if there is single rRNA FASTA file present in this folder

if [[ ${TEMP2} -gt 0 ]]; then
  log "Total of ${TEMP2} rRNA FASTA files found in '${R_RNA_PATH}'"
  log "Skipping copying FATSA files from SortMeRNA database"
else
  log "No rRNA FATSA files found. Trying to copy from SortMeRNA database"

  # Assuming all rRNA database is in 'rRNA_databases' folder
  find "${TOOL_SORTMERNA}/rRNA_databases/" -name "*.fasta" |
    while read -r f; do
      cp "${f}" "${R_RNA_PATH}" # Copy to our reference folder
    done

  log "rRNA FATSA files copied from SortMeRNA database"
fi

log "Starting rRNA indexing"

find "${R_RNA_PATH}" -maxdepth 1 -name "*.fasta" |
  while read -r f; do
    name=$(echo "${f}" | xargs -I {} basename {}) # Get file name
    # Check if indexing has done already
    if [[ -n $(find "${R_RNA_PATH}/index/" -name "${name}*") ]]; then
      log "Indexing has already done for ${name}"
    else
      index_rrna "${f}" "${name}"
    fi
  done
