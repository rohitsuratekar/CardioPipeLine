#!/usr/bin/env bash
#===============================================================================
#
#          FILE: rrna_filtering.sh
#
#         USAGE: ./rrna_filtering.sh
#
#   DESCRIPTION: Filters the sequence by seperating rRNA from regular RNA
#
#  REQUIREMENTS: rRNA indexing
#        AUTHOR: Rohit Suratekar
#  ORGANIZATION: IIMCB
#       CREATED: Friday 09 August 2019 16:26
#      REVISION: 3
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

# Create local merge folder
merge_folder="${FOLDER_FILTERING}/merged"
mkdir -p "$merge_folder"

# This analysis assume Paired-end SRA file which has split into 2 parts

readonly FORWARD_READS="${FOLDER_DATA}/${SRA_ID}_1.fastq"
readonly REVERSE_READS="${FOLDER_DATA}/${SRA_ID}_2.fastq"

perform_merge() {
  # Merges two files to create one as required for SortMeRNA

  "${TOOL_SORTMERNA}"/scripts/merge-paired-reads.sh \
    "${FORWARD_READS}" "${REVERSE_READS}" \
    "${1}/${SRA_ID}${READ_MERGED_EXTENSION}"
}

# To avoid ambiguty, name of this function is kept different
unmerge_sequence() {
  # Unmerges the sequences for further processing
  "${TOOL_SORTMERNA}"/scripts/unmerge-paired-reads.sh \
    "${FOLDER_FILTERING}/${SRA_ID}${EXTENSION_R_RNA_FILTERED}.fastq" \
    "${FOLDER_FILTERING}/${SRA_ID}${EXTENSION_R_RNA_FILTERED}_1.fastq" \
    "${FOLDER_FILTERING}/${SRA_ID}${EXTENSION_R_RNA_FILTERED}_2.fastq"
}

# Check if merged folder already exists

if check_file -p "$merge_folder" "${SRA_ID}${EXTENSION_MERGED_READS}"; then
  log "Merged file already exists in '${merge_folder}'"
else
  if ! perform_merge "$merge_folder"; then
    log "Unable to merge the fastq files for rRNA filtering"
  fi
fi

# Make folder to keep aligned sequences
mkdir -p "${FOLDER_FILTERING}"
# Make Folder for Key-value data-store

KVDS_PATH="${FOLDER_FILTERING}/kvds"
rm -rf "${KVDS_PATH}" # This is needed because sortmerna script won't run if this folder is non-empty
mkdir -p "${KVDS_PATH}"

filter_rrna() {

  # Start SortMeRNA
  #
  # --ref : rRNA reference sequences and their indices
  # --reads : Merged Fasta file
  # --aligned : Aligned file name and path
  # --other : non-aligned file name and path (this is non-rRNA)
  # --fastx : fasta format
  # --log : Creates log files
  # -v : Verbose output
  # -a : number of threads (if you don't know how many core your system has, you can ommit this option)
  # -d : key value data store dictionary path (This option is needed because sortmerna won't delete kvds folder file by itself)
  # --paired_in : If either one from forward or reverse read matches, it will be in aligned file

  "${TOOL_SORTMERNA}"/bin/sortmerna \
    --ref "${2}" \
    --reads "${3}/${SRA_ID}${EXTENSION_MERGED_READS}" \
    --aligned "${1}/${SRA_ID}${EXTENSION_R_RNA}" \
    --other "${1}/${SRA_ID}${EXTENSION_R_RNA_FILTERED}" \
    --fastx \
    --paired_in \
    --log \
    -a 8 \
    -v \
    -d "${KVDS_PATH}"

}

# Generate sequence that can be used with SortMeRNA command

REF_SEQ=""
while read -r f; do
  name=$(echo "${f}" | xargs -I {} basename {}) # Get file name
  REF_SEQ="${REF_SEQ}${f},${FOLDER_R_RNA}/index/${name}${EXTENSION_MERGED_READS}:"
done <<<"$(find "${FOLDER_R_RNA}" -maxdepth 1 -name "*.fasta" -type f -size +1M)"

# Remove last character
REF_SEQ="echo ${REF_SEQ/.$/}"

# Skip filtering if it is already done
if check_file -p "$FOLDER_FILTERING" "${SRA_ID}${EXTENSION_R_RNA_FILTERED}.fastq"; then
  log "Filtered file already exist. Skipping rRNA filtering"
else
  log "starting rRNA filtering"
  filter_rrna "$FOLDER_FILTERING" "$REF_SEQ" "$merge_folder"
fi

# Check if unmerging has done or not
# if not, perform the unmerging


if check_file -p "$FOLDER_FILTERING" "${SRA_ID}${EXTENSION_R_RNA_FILTERED}_*.fastq"; then
  log "Filtered-split sequences already exists in ${FOLDER_FILTERING}"
  log "Skipping splitting of filtered sequence"
else
  log "Splitting filtered sequence"
  if unmerge_sequence ; then
    log "Splitting completed in folder ${FOLDER_FILTERING}"
  else
    log "Problem encountered while unmerging the sequence"
    exist 1
  fi
fi
