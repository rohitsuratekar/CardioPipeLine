import os

BASE = config['base']


def get_fastq(wildcards):
    base = f"{BASE}/fastq/{wildcards.SRR_ID}"
    if layout_expand(wildcards.SRR_ID):
        return [f"{base}/{wildcards.SRR_ID}.sra_1.fastq",
                f"{base}/{wildcards.SRR_ID}.sra_2.fastq"]
    else:
        return [f"{base}/{wildcards.SRR_ID}.sra.fastq"]


rule convert_fastq_paired:
    input:
         fasterq=config["tools"]["fasterq-dump"],
         sra_file=expand("{folder}/srr/{srr}.sra", folder=BASE,
                         srr="{SRR_ID}")
    output:
          expand("{folder}/fastq/{srr}/{srr}.sra_{i}.fastq", folder=BASE,
                 srr="{SRR_ID}", i=[1, 2])
    threads: config['threads']
    # Remember to keep space after each of the options
    shell:
         "{input.fasterq} "  # Binary
         "{input.sra_file} "  # Input SRA File
         "--split-files "  # Split files according to the reads
         "--progres "  # Show progress
         "--skip-technical "  # Skip technical reads
         "--threads {threads} "  # Number of threads
         "--outdir {BASE}/fastq/{wildcards.SRR_ID}"

rule convert_fastq_single:
    input:
         fasterq=config["tools"]["fasterq-dump"],
         sra_file=expand("{folder}/srr/{srr}.sra", folder=BASE,
                         srr="{SRR_ID}")
    output:
          expand("{folder}/fastq/{srr}/{srr}.sra.fastq", folder=BASE,
                 srr="{SRR_ID}")
    threads: config['threads']
    # Remember to keep space after each of the options
    shell:
         "{input.fasterq} "  # Binary
         "{input.sra_file} "  # Input SRA File
         "--split-files "  # Split files according to the reads
         "--progres "  # Show progress
         "--skip-technical "  # Skip technical reads
         "--threads {threads} "  # Number of threads
         "--outdir {BASE}/fastq/{wildcards.SRR_ID}"


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

rule filter_rrna_paired:
    input:
         reads=get_fastq,
         refs=get_refs,
         sortmerna=config["tools"]["sortmerna"]
    output:
          expand("{folder}/filtered/{srr}/{srr}.filtered_{i}.fastq",
                 folder=BASE, srr="{SRR_ID}", i=[1, 2])
    log:
       expand("{folder}/logs/{srr}.filtered.log", folder=BASE, srr="{SRR_ID}")
    threads: config['threads']
    params:
          ref_str=lambda wildcards, input: " ".join(
              ["-ref {}".format(x) for x in input.refs]),
          read_str=lambda wildcards, input: " ".join(
              ["-reads {}".format(x) for x in input.reads]),
          index_dir="{}/index/sortmerna".format(BASE),
          filter_pre="{SRR_ID}.filtered_",
    shell:
         """
         rm -rf {params.index_dir}/kvdb
         
         {input.sortmerna} {params.ref_str} {params.read_str} -workdir {params.index_dir} -fastx -other -a {threads} -num_alignments 1 -paired_in
         
         mkdir -p {BASE}/filtered/{wildcards.SRR_ID}
         
         file1="{BASE}/filtered/{wildcards.SRR_ID}/{params.filter_pre}1.fastq"
         file2="{BASE}/filtered/{wildcards.SRR_ID}/{params.filter_pre}2.fastq"  
         scripts/unmerge-paired-reads.sh {params.index_dir}/out/other.fastq $file1 $file2 
         
         mv {params.index_dir}/out/aligned.log {log} 
         rm -f {params.index_dir}/out/aligned.fastq
         rm -f {params.index_dir}/out/other.fastq
         """

rule filter_rrna_single:
    input:
         reads=get_fastq,
         refs=get_refs,
         sortmerna=config["tools"]["sortmerna"]
    output:
          expand("{folder}/filtered/{srr}/{srr}.filtered.fastq",
                 folder=BASE, srr="{SRR_ID}")
    log:
       expand("{folder}/logs/{srr}.filtered.log", folder=BASE, srr="{SRR_ID}")
    threads: config['threads']
    params:
          ref_str=lambda wildcards, input: " ".join(
              ["-ref {}".format(x) for x in input.refs]),
          read_str=lambda wildcards, input: " ".join(
              ["-reads {}".format(x) for x in input.reads]),
          index_dir="{}/index/sortmerna".format(BASE),
    shell:
         """
         rm -rf {params.index_dir}/kvdb
         
         {input.sortmerna} {params.ref_str} {params.read_str} -workdir {params.index_dir} -fastx -other -a {threads} -num_alignments 1 -paired_in
         
         mkdir -p {BASE}/filtered/{wildcards.SRR_ID} 
         mv {params.index_dir}/out/other.fastq {output}
         
         mv {params.index_dir}/out/aligned.log {log} 
         rm -f {params.index_dir}/out/aligned.fastq
         rm -f {params.index_dir}/out/other.fastq
         """


def method_reads(wildcards):
    if config["filter_rrna"] == 1:
        base = f"{BASE}/filtered/{wildcards.SRR_ID}"
        if layout_expand(wildcards.SRR_ID):
            return [f"{base}/{wildcards.SRR_ID}.filtered_1.fastq",
                    f"{base}/{wildcards.SRR_ID}.filtered_2.fastq"]
        else:
            return [f"{base}/{wildcards.SRR_ID}.filtered.fastq"]

    else:
        return get_fastq(wildcards)
