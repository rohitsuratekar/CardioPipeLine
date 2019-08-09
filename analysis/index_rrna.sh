#!/usr/bin/env bash
#===============================================================================
#
#          FILE: index_rrna.sh
#
#         USAGE: ./index_rrna.sh
#
#   DESCRIPTION: Indexes rRNA dataset which can be used for SortMeRNA
#
#  REQUIREMENTS: rRNA fasta files
#         NOTES: Usually rRNA files are kept in SortMeRNA's 'rRNA_databases'
#                folder. You can use those default files or can use your own by
#                copying them to $FOLDER_R_RNA
#        AUTHOR: Rohit Suratekar
#  ORGANIZATION: IIMCB
#       CREATED: Friday 09 August 2019 15:58
#      REVISION:  3
#===============================================================================

set -o nounset # Treat unset variables as an error

source config.sh

# If it is configure to skip the rRNA filtering, then just retuen
# the script as a successful

if [[ ${RNA_FILTERING} -eq 0 ]]; then
  log "Skipping rRNA indexing. If you want to filter sequences with rRNA
  then change RNA_FILTERING value to 1 in config.sh file."
  exit 0
fi

# Index files will be stored here
mkdir -p "$FOLDER_R_RNA_INDEX"

index_rrna() {
  # Following command is VERY sensitive to the spaces. So double check if
  # you are not giving any additional spaces in the command

  "${TOOL_SORTMERNA}/bin/indexdb" --ref "${1}","${FOLDER_R_RNA}/index/${2}${EXTENSION_R_RNA_INDEX}" -v # Make Index

}

copy_rrna_files() {

  # Check if any of the fasta files are avilable in the default location
  if ! check_file -p "${TOOL_SORTMERNA}/rRNA_databases/" "*.fasta"; then
    log "No rRNA file found in '${TOOL_SORTMERNA}/rRNA_databases/'"
    log "For rRNA filtering, these files are important."
    log "You can manually copy rRNA fasta files to ${FOLDER_R_RNA}"
    exit 1
  fi

  # Copies rRNA files from SortMeRNA tool folder to desired folder
  find "${TOOL_SORTMERNA}/rRNA_databases/" -name "*.fasta" |
    while read -r f; do
      cp "${f}" "${FOLDER_R_RNA}" # Copy to our reference folder
    done
  log "rRNA FATSA files copied from SortMeRNA database"
}

# Check if rRNA FASTA files exists
if check_file -p "$FOLDER_R_RNA" "*.fasta"; then
  log "rRNA FASTA files found in '${FOLDER_R_RNA}'"
  log "Skipping copying FATSA files from SortMeRNA database"
else
  log "No rRNA FATSA files found. Trying to copy from SortMeRNA database"
  copy_rrna_files
fi

# Start the rRNA indexing

find "${FOLDER_R_RNA}" -maxdepth 1 -name "*.fasta" |
  while read -r f; do
    name=$(echo "${f}" | xargs -I {} basename {}) # Get file name
    # Check if indexing has done already
    if [[ -n $(find "${FOLDER_R_RNA}/index/" -name "${name}*") ]]; then
      log "Indexing has already done for ${name}"
    else
      index_rrna "${f}" "${name}"
    fi
  done
