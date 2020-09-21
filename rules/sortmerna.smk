"""
CardioPipeLine (c) 2020

Author: Rohit Suratekar, ZDG Lab, IIMCB

This deals with all rules related SortMeRNA tools and rRNA-filtering
"""

import os

BASE = config['base']


def get_refs(wildcards):
    path = config['genome']['rrna_db']
    files = []
    for f in os.listdir(path):
        if f.endswith("fasta") and f.startswith("silva-euk"):  # Use only
            # eukaryotic data
            files.append(f"{path}/{f}")
    return files


# Option explanation
# All references for indexing will be input as:  -ref ref1.fa -ref ref2.fa ...
# All samples will be input as: -reads read1.fastq -reads read2.fastq
# -workdir : Where index will be kept
# -fastx : Outputs only in fastq format
# -other : Output the Other.fastq file (which has NON-rRNA reads)
# -a : Number of threads
# -num_alignment :  Only first alignment passing of E value (fastest option
# recommended in the manual is to set this value to 1)
# -paired_in : if either of read match, put both reads into the 'aligned.fasta' file
# --out2 : Fastq files will be split automatically

# For Paired layout
rule filter_rrna_paired:
    input:
         reads=["{BASE_FOLDER}/fastq/{SRR_ID}.sra_1.fastq",
                "{BASE_FOLDER}/fastq/{SRR_ID}.sra_2.fastq"],
         refs=get_refs
    output:
          ["{BASE_FOLDER}/filtered/{SRR_ID}.sra.filtered_1.fastq",
           "{BASE_FOLDER}/filtered/{SRR_ID}.sra.filtered_2.fastq"]

    threads: config['threads']
    params:
          sortmerna=config["tools"]["sortmerna"],
          ref_str=lambda wildcards, input: " ".join(
              ["-ref {}".format(x) for x in input.refs]),
          read_str=lambda wildcards, input: " ".join(
              ["-reads {}".format(x) for x in input.reads]),
          index_dir="{}/index/sortmerna".format(BASE),
    shell:
         """
         rm -rf {params.index_dir}/kvdb
         
         {params.sortmerna} {params.ref_str} {params.read_str} -workdir {params.index_dir} -fastx -other -a {threads} -num_alignments 1 -paired_in --out2 
         
         mkdir -p {BASE}/filtered/ 
         
         mkdir -p {BASE}/logs
         
         mv {params.index_dir}/out/aligned.log {BASE}/logs/{wildcards.SRR_ID}_filtering.log
         
         mv {params.index_dir}/out/other_fwd.fastq {BASE}/filtered/{wildcards.SRR_ID}.sra.filtered_1.fastq
         mv {params.index_dir}/out/other_rev.fastq {BASE}/filtered/{wildcards.SRR_ID}.sra.filtered_2.fastq
         
         rm  -f {params.index_dir}/out/aligned_fwd.fastq
         rm  -f {params.index_dir}/out/aligned_rev.fastq
         """

# For single layout
rule filter_rrna_single:
    input:
         reads="{BASE_FOLDER}/fastq/{SRR_ID}.sra.fastq",
         refs=get_refs
    output:
          "{BASE_FOLDER}/filtered/{SRR_ID}.sra.filtered.fastq"
    threads: config['threads']
    params:
          sortmerna=config["tools"]["sortmerna"],
          ref_str=lambda wildcards, input: " ".join(
              ["-ref {}".format(x) for x in input.refs]),
          index_dir="{}/index/sortmerna".format(BASE),
    shell:
         """
         rm -rf {params.index_dir}/kvdb
         
         {params.sortmerna} {params.ref_str}  -reads {input.reads} -workdir {params.index_dir} -fastx -other -a {threads} -num_alignments 1 -paired_in 
         
         mkdir -p {BASE}/filtered/ 
         
         mkdir -p {BASE}/logs
         
         mv {params.index_dir}/out/aligned.log {BASE}/logs/{wildcards.SRR_ID}_filtering.log
         
         mv {params.index_dir}/out/other.fastq {BASE}/filtered/{wildcards.SRR_ID}.sra.filtered.fastq 
         
         rm  -f {params.index_dir}/out/aligned.fastq
         """
