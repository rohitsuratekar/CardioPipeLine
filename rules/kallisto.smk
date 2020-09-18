"""
CardioPipeLine (c) 2020

Author: Rohit Suratekar, ZDG Lab, IIMCB

This deals with all rules related to Kallisto
"""

rule index_kallisto:
    input:
         kallisto=config["tools"]["kallisto"],
         transcript=ancient(config["genome"]["transcript"])

    output: "{BASE}/index/kallisto/kallisto.idx"

    shell:
         "{input.kallisto} index "  # Kallisto index mode
         "-i {output} "  # Filename (not folder) for kallisto
         "{input.transcript}"  # Transcript


def make_kallisto_reads(files):
    if len(files) == 2:
        # Paired end reads
        return " ".join(files)
    else:
        # Single end, use parameters from config file
        frg_length = config["extra"]["kallisto-fragment-length"]  # Fragment
        # length
        frg_sd = config["extra"]["kallisto-standard-deviation"]  # Standard
        # deviation
        return f"--single -l {frg_length} -s {frg_sd} {files}"


rule kallisto_quant:
    input:
         kallisto=config["tools"]["kallisto"],
         k_index=ancient("{BASE}/index/kallisto/kallisto.idx"),
         gtf=ancient(config["genome"]["gtf_annotation"]),
         files=get_fastq_files,
    threads: config["threads"]
    output:
          "{BASE}/mappings/kallisto/{SRR_ID}/abundance.tsv"
    params:
          reads=lambda wildcards, input: make_kallisto_reads(input.files)
    shell:
         "{input.kallisto} quant "  # Kallisto quant mode
         "-i {input.k_index} "  # Path to the index file
         "-o {wildcards.BASE}/mappings/kallisto/{wildcards.SRR_ID} "  # Output 
         # folder
         "-g {input.gtf} "  # Genome GTF file
         "-t {threads} "  # Number of threads
         "--bias "  # Sequence based bias correction,
         "--plaintext "  # Plaintext output instead HDF5
         "{params.reads}"  # Add reads and layout specific options
