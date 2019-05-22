#!/usr/bin/env bash
#Creates index files from given rRNA fasta files

echo Script Started at : $(date)

# Path of folder where you want your fasta files and index files to be stored
readonly rna_path="/home/rsuratekar/Documents/Rohit/Genome/rrna/"

# Index files will be stored here
mkdir ${rna_path}index

# Path of SortMeRNA folder
readonly tool_path="/home/rsuratekar/Documents/Rohit/Tools/sortmerna-3.0.3/"

# Assuming all rRNA database is in 'rRNA_databases' folder
find "${tool_path}rRNA_databases/" -name "*.fasta"|
while read f; 
do	
	cp $f $rna_path
	name=$(echo $f| xargs -I {} basename {})
	${tool_path}bin/indexdb --ref $f,${rna_path}index/${name}.idx -v
done
