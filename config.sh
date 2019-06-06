#!/usr/bin/env bash
# All variables related to RNA-seq pipeline can be changed here.

source ./local.sh  # For importing local variable. You can remove this line and use your own "ORIGIN_LOCAL"


SRA_ID="SRX4720625" #ID of SRA file. This will be used in many output and input names
GENOME_FASTA="Danio_rerio.GRCz11.dna.primary_assembly.fa" #Name of the fasta reference file
GENOME_ANNOTATION="Danio_rerio.GRCz11.96.chr.gtf" # Genome Annotation file
R_RNA_INDEX_EXTENSION=".idx" # Extension which will be added at the end of rRNA index file
R_RNA_FILTERED_EXTENSION=".filtered" # Extension which will be added at the end of aligned rRNA filtered sequence
READ_MERGED_EXTENSION="_merged.fastq" # Extension will be added at the end of merged file

# Log related options
LOG="log.txt" # Name of the log file
LOG_ENABLED=0 # 1=Yes and 0=NO (Add to Log File)
LOG_DISPLAY=1 # 1=Yes and 0=No (Display on Terminal)

# Various paths
ORIGIN=${ORIGIN_LOCAL} # Base folder for everything e.g /home/user/path/to/base/folder
GENOME_HOME="$ORIGIN/Genome" # Standard genomes and their indexes
TOOL="$ORIGIN/Tools" # Tools involved in this pipeline
PIPE="$ORIGIN/Pipeline" # All pipeline scripts
DATA="$ORIGIN/Data" # All experimental data files
NCBI_DOWNLOAD=${NCBI_LOCAL} # /home/user/ncbi/public/sra
RAW_DATA="$DATA/raw" # Here all downloaded raw data will be stored

mkdir -p ${RAW_DATA} # Make folder for raw data if not present


# All Tools paths

TOOL_SRA="$TOOL/sratoolkit.2.9.6-1-ubuntu64"
TOOL_HISAT="$TOOL/hisat2-2.1.0"
TOOL_HT_SEQ="$TOOL/htseq-release_0.11.0"
TOOL_SORTMERNA="$TOOL/sortmerna-3.0.3"
TOOL_FASTQC="$TOOL/FastQC"
TOOL_STAR="$TOOL/STAR-2.7.1a/source"
TOOL_SAMTOOL="$TOOL/samtools-1.9"
TOOL_STRINGTIE="$TOOL/stringtie-1.3.4d.Linux_x86_64"


# Log Function

log(){
# Log the details of it is enabled
if [[ ${LOG_ENABLED} -eq 1 ]]
then
	echo "RLOG $(date +%y-%m-%d" "%H:%M:%S) :$1" >> ${LOG}
fi
# Display if display is enabled
if [[ ${LOG_DISPLAY} -eq 1 ]]
then
	echo "RLOG $(date +%y-%m-%d" "%H:%M:%S) :$1"
fi
}


log "===================================="
log "Script: $(basename $0)"
log "Time: $(date)"
log "===================================="