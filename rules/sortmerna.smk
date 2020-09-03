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

rule filter_rrna_paired:
    input:
         reads=expand("{folder}/fastq/{srr}.sra_{i}.fastq",
                      folder=config['base'], srr="{SRR_ID}", i=[1, 2]),
         refs=get_refs
    output:
          expand("{folder}/filtered/{srr}.sra.filtered_{i}.fastq",
                 folder=config['base'], srr="{SRR_ID}", i=[1, 2])

    threads: config['threads']
    params:
          sortmerna=config["tools"]["sortmerna"],
          ref_str=lambda wildcards, input: " ".join(
              ["-ref {}".format(x) for x in input.refs]),
          read_str=lambda wildcards, input: " ".join(
              ["-reads {}".format(x) for x in input.reads]),
          index_dir="{}/index/sortmerna".format(BASE),
          filter_pre=lambda wildcards: "{}.filtered_".format(wildcards.SRR_ID),
    shell:
         """
         rm -rf {params.index_dir}/kvdb
         
         {params.sortmerna} {params.ref_str} {params.read_str} -workdir {params.index_dir} -fastx -other -a {threads} -num_alignments 1 -paired_in --out2 
         
         mkdir -p {BASE}/filtered/ 
         
         mkdir -p {BASE}/logs
         
         mv {params.index_dir}/out/aligned.log {BASE}/logs/{wildcards.SRR_ID}_filtering.log
         """
