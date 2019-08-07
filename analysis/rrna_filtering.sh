#!/usr/bin/env bash
#
# Filter sequence with rRNA contamination

source config.sh

# If it is configure to skip the rRNA filtering, then just retuen
# the script as a successful

if [[ ${RNA_FILTERING} -eq 0 ]]; then
  log "Skipping rRNA indexing. If you want to filter sequences with rRNA
  then change RNA_FILTERING value to 1 in config.sh file."
  exit 0
fi

log "Checking if merged file required for filtering exists"

readonly MERGE_FOLDER="${DATA}/filtering/merged"
readonly R_RNA_PATH="$GENOME_HOME/rrna"

mkdir -p "${MERGE_FOLDER}"

# This analysis assume Paired-end SRA file which has split into 2 parts

readonly FORWARD_READS="${DATA}/${SRA_ID}_1.fastq"
readonly REVERSE_READS="${DATA}/${SRA_ID}_2.fastq"

if [[ -n $(find "${MERGE_FOLDER}" -name "${SRA_ID}${READ_MERGED_EXTENSION}" -size +1M) ]]; then
  log "Merged file already exists in '${MERGE_FOLDER}'"
else
  "${TOOL_SORTMERNA}"/scripts/merge-paired-reads.sh \
    "${FORWARD_READS}" "${REVERSE_READS}" \
    "${MERGE_FOLDER}/${SRA_ID}${READ_MERGED_EXTENSION}"
  # Check if merging is successful
  if [[ "$?" -ne 0 ]]; then
    log "Unable to merge fastq files for rRNA filtering. Exiting script"
    exit 1
  else
    log "Merged reads generated for rRNA filtering"
  fi
fi

# Generate sequence that can be used with SortMeRNA command

REF_SEQ=""
while read -r f; do
  name=$(echo "${f}" | xargs -I {} basename {}) # Get file name
  REF_SEQ="${REF_SEQ}${f},${R_RNA_PATH}/index/${name}${R_RNA_INDEX_EXTENSION}:"
done <<<"$(find "${R_RNA_PATH}" -maxdepth 1 -name "*.fasta" -type f -size +1M)"

# Remove last character
REF_SEQ="$(echo "${REF_SEQ/.$/}")"

log "rRNA filtering initiating"

# Make folder to keep aligned sequences

mkdir -p "${FOLDER_FILTERED}"

# Make Folder for Key-value data-store

KVDS_PATH="${DATA}/filtering/kvds"
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

  if [[ -n $(find "${1}" -maxdepth 1 -type f -name "${SRA_ID}${R_RNA_FILTERED_EXTENSION}.fastq" -size +0M) ]]; then
    log "Filtered file already exist. Skipping rRNA filtering"
  else
    log "starting rRNA filtering"

    "${TOOL_SORTMERNA}"/bin/sortmerna \
      --ref "${2}" \
      --reads "${3}/${SRA_ID}${READ_MERGED_EXTENSION}" \
      --aligned "${1}/${SRA_ID}_rRNA" \
      --other "${1}/${SRA_ID}${R_RNA_FILTERED_EXTENSION}" \
      --fastx \
      --paired_in \
      --log \
      -a 8 \
      -v \
      -d "${KVDS_PATH}"

    # Check if command is successful
    if [[ "$?" -ne 0 ]]; then
      log "Unable to filter sequence with rRNA. Exiting script"
      exit 1
    else
      log "rRNA filtering completed"
    fi

  fi

}

# Start filtering
filter_rrna "${FOLDER_FILTERED}" "${REF_SEQ}" "${MERGE_FOLDER}"

# Split the merged file from SortMeRNA
TEMP3=$(find "${FOLDER_FILTERED}" -maxdepth 1 -type f -name "${SRA_ID}${R_RNA_FILTERED_EXTENSION}_*.fastq" -size +0M | wc -l)

if [[ ${TEMP3} -eq 2 ]]; then
  log "${TEMP3} filtered-split sequences already exists in ${FOLDER_FILTERED}"
  log "Skipping splitting of filtered sequence"
else
  log "Splitting filtered sequence"

  "${TOOL_SORTMERNA}"/scripts/unmerge-paired-reads.sh \
    "${FOLDER_FILTERED}/${SRA_ID}${R_RNA_FILTERED_EXTENSION}.fastq" \
    "${FOLDER_FILTERED}/${SRA_ID}${R_RNA_FILTERED_EXTENSION}_1.fastq" \
    "${FOLDER_FILTERED}/${SRA_ID}${R_RNA_FILTERED_EXTENSION}_2.fastq"

  # Check if command is successful
  if [[ "$?" -ne 0 ]]; then
    log "Unable to filter sequence with rRNA. Exiting script"
    exit 1
  else
    log "Splitting completed in folder ${FOLDER_FILTERED}"
  fi

fi
