"""
CardioPipeLine (c) 2020

Author: Rohit Suratekar, ZDG Lab, IIMCB

This deals with all rules related to FastQC
"""

# -t : Threads
# -o : OutDir

rule quality_paired:
    input:
         f1="{BASE_FOLDER}/fastq/{SRR_ID}.sra_1.fastq",
         f2="{BASE_FOLDER}/fastq/{SRR_ID}.sra_2.fastq",
         fastqc=config["tools"]["fastqc"]
    threads: config["threads"]
    output:
          o1="{BASE_FOLDER}/quality/{SRR_ID}/{SRR_ID}.sra_1_fastqc.html",
          o2="{BASE_FOLDER}/quality/{SRR_ID}/{SRR_ID}.sra_2_fastqc.html"
    shell:
         """
         {input.fastqc} -t {threads} -o {wildcards.BASE_FOLDER}/quality/{wildcards.SRR_ID} {input.f1} {input.f2}
         """

rule quality_single:
    input:
         f1="{BASE_FOLDER}/fastq/{SRR_ID}.sra.fastq",
         fastqc=config["tools"]["fastqc"]
    threads: config["threads"]
    output:
          o1="{BASE_FOLDER}/quality/{SRR_ID}/{SRR_ID}.sra_fastqc.html"
    shell:
         """
         {input.fastqc} -t {threads} -o {wildcards.BASE_FOLDER}/quality/{wildcards.SRR_ID} {input.f1}
         """

rule quality_filtered_paired:
    input:
         f1="{BASE_FOLDER}/filtered/{SRR_ID}.sra.filtered_1.fastq",
         f2="{BASE_FOLDER}/filtered/{SRR_ID}.sra.filtered_2.fastq",
         fastqc=config["tools"]["fastqc"]
    threads: config["threads"]
    output:
          o1="{BASE_FOLDER}/quality/{SRR_ID}/{SRR_ID}.sra.filtered_1_fastqc.html",
          o2="{BASE_FOLDER}/quality/{SRR_ID}/{SRR_ID}.sra.filtered_2_fastqc.html",
    shell:
         """
         {input.fastqc} -t {threads} -o {wildcards.BASE_FOLDER}/quality/{wildcards.SRR_ID} {input.f1} {input.f2}
         """

rule quality_filtered_single:
    input:
         f1="{BASE_FOLDER}/filtered/{SRR_ID}.sra.filtered.fastq",
         fastqc=config["tools"]["fastqc"]
    threads: config["threads"]
    output:
          o1="{BASE_FOLDER}/quality/{SRR_ID}/{SRR_ID}.sra.filtered_fastqc.html"
    shell:
         """
         {input.fastqc} -t {threads} -o {wildcards.BASE_FOLDER}/quality/{wildcards.SRR_ID} {input.f1}
         """
