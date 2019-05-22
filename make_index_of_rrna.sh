#!/usr/bin/env bash
#Creates index files from given rRNA fasta files

source ./config.sh

# Path of folder where you want your fasta files and index files to be stored
readonly rna_path="$GENOME_HOME/rrna/"

# Index files will be stored here
mkdir ${rna_path}index

# Path of SortMeRNA folder
readonly tool_path="$TOOL_SORTMERNA/"

# Assuming all rRNA database is in 'rRNA_databases' folder
find "${tool_path}rRNA_databases/" -name "*.fasta"|
while read f; 
do	
	cp ${f} ${rna_path}
	name=$(echo ${f}| xargs -I  basename )
	${tool_path}bin/indexdb --ref ${f},${rna_path}index/${name}.idx -v
done
