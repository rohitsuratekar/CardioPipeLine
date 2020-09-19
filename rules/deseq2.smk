"""
CardioPipeLine (c) 2020

Author: Rohit Suratekar, ZDG Lab, IIMCB

This deals with all rules related to DESeq2 analysis
"""


# Important: Script rules are relative to this file instead working dir.
# Hence we need to use '../scripts'

def generate_all_bams(wildcards):
    files = []
    for _, row in SAMPLES_DF.iterrows():
        srr = row["run"]
        files.append(
            f"{config['base']}/bams/{srr}.sra.Aligned.sortedByCoord.out.bam")
    return files


rule generate_counts:
    input:
         samples=config["samples"],
         gtf=config["genome"]["gtf_annotation"],
         files=generate_all_bams

    output: "{BASE}/test.txt"
    params:
          mode="star"
    script: "../scripts/deseq2.R"
