"""
CardioPipeLine (c) 2020

Author: Rohit Suratekar, ZDG Lab, IIMCB

This deals with all rules related to Salmon
"""

# Salmon index will be generated with the decoy file.
# First generate that decoy file
rule generate_decoy:
    input:
         gtf=ancient(config["genome"]["gtf_annotation"]),
         fasta=ancient(config["genome"]["fasta"]),
         transcript=ancient(config["genome"]["transcript"]),
         bedtools=config["tools"]["bedtools"],
         mashmap=config["tools"]["mashmap"]
    output:
          out_index="{BASE}/index/salmon/gentrome.fa",
          out_keys="{BASE}/index/salmon/decoys.txt"
    threads: config["threads"]
    shell:
         "scripts/generateDecoyTranscriptome.sh "
         "-j {threads} "  # Number of threads
         "-a {input.gtf} "  # Genome GTF File
         "-g {input.fasta} "  # Genome in fasta
         "-t {input.transcript} "  # Transcript file in fasta
         "-o {wildcards.BASE}/index/salmon "  # Salmon index output folder
         "-b {input.bedtools} "  # Path to bedtools
         "-m {input.mashmap} "  # Path to mashmap

# Generate index for salmon
rule index_salmon:
    input:
         decoy_index=ancient("{BASE}/index/salmon/gentrome.fa"),
         decoy_keys=ancient("{BASE}/index/salmon/decoys.txt"),
         salmon=config["tools"]["salmon"]
    threads: config["threads"]
    output: "{BASE}/index/salmon/pos.bin"
    shell:
         "{input.salmon} index "  # Salmon index mode
         "-t {input.decoy_index} "  # Decoy transcript
         "-i {wildcards.BASE}/index/salmon "  # Path to index folder
         "--decoys {input.decoy_keys} "  # Decoy keys
         "-k 31 "  # K-mer (31) is recommended by the authors
         "-p {threads}"  # Number of threads


def get_salmon_reads(files):
    if isinstance(files, str):
        return f"-r {files}"
    else:
        sorted_list = sorted(files)
        return f"-1 {sorted_list[0]} -2 {sorted_list[1]}"


# Perform Salmon mapping
rule salmon_quant:
    input:
         salmon=config["tools"]["salmon"],
         salmon_index=ancient("{BASE}/index/salmon/pos.bin"),  # This is
         # just needed to ensure that index is prepared.
         gtf=config["genome"]["gtf_annotation"],
         files=get_fastq_files,
    output:
          "{BASE}/mappings/salmon/{SRR_ID}/quant.sf"

    params:
          reads=lambda wildcards, input: get_salmon_reads(input.files)

    threads: config["threads"]
    shell:
         "{input.salmon} quant "  # Salmon quant mode
         "-i {wildcards.BASE}/index/salmon "  # Path of the index
         "-l A "  # Automatic library type
         "--validateMapping "  #enables selective alignment of the 
         # sequencing reads
         "--gcBias "  # Corrects for  fragment-level GC biases. Also 
         # recommended by the DESEq2 analysis
         "-p {threads} "  # Number of threads   
         "-g {input.gtf} "  # Additional GTF file to get quants.genes.sf
         "-o {wildcards.BASE}/mappings/salmon/{wildcards.SRR_ID} "  # output 
         # folder
         "{params.reads}"  # Input reads
