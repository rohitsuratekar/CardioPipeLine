import os
import shutil

BASE = config['base']

checkpoint convert_fastq:
    input:
         fasterq=config["tools"]["fasterq-dump"],
         sra_file=expand("{folder}/srr/{srr}.sra", folder=BASE,
                         srr="{SRR_ID}")
    output: out=directory(os.path.join(BASE, "fastq", "{SRR_ID}"))
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


def get_fastq(wildcards):
    base = "{}/fastq/{}".format(BASE, wildcards.SRR_ID)
    if os.path.exists(base):
        files = glob_wildcards(os.path.join(base, "{name}.fastq")).name
        if len(files) == 0:
            # Folder is empty and remove it
            shutil.rmtree(base)
            print("Deleted empty folder for fastq conversion")

    path = checkpoints.convert_fastq.get(**wildcards).output.out
    return expand("{path}/{file}.fastq",
                  path=path,
                  file=glob_wildcards(os.path.join(path, "{name}.fastq")).name)


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

rule filter_rrna:
    input:
         reads=get_fastq,
         refs=get_refs,
         sortmerna=config["tools"]["sortmerna"]
    output: directory(BASE + "/filtered/{SRR_ID}")
    log:
       expand("{folder}/logs/{srr}.filtered.log", folder=BASE, srr="{SRR_ID}")
    threads: config['threads']
    params:
          ref_str=lambda wildcards, input: " ".join(
              ["-ref {}".format(x) for x in input.refs]),
          read_str=lambda wildcards, input: " ".join(
              ["-reads {}".format(x) for x in input.reads]),
          index_dir="{}/index/sortmerna".format(BASE),
          paired=lambda wildcards, input: int(len(input.reads) == 2),
          filter_pre="{SRR_ID}.filtered_",
    shell:
         """
         rm -rf {params.index_dir}/kvdb
         
         {input.sortmerna} {params.ref_str} {params.read_str} -workdir {params.index_dir} -fastx -other -a {threads} -num_alignments 1 -paired_in
         
         mkdir -p {output}
         
         if [[ {params.paired} -eq 1 ]]; then
            file1="{output}/{params.filter_pre}1.fastq"
            file2="{output}/{params.filter_pre}2.fastq"  
            scripts/unmerge-paired-reads.sh {params.index_dir}/out/other.fastq $file1 $file2 
         else
            file3="{output}/{wildcards.SRR_ID}.filtered.fastq" 
            mv {params.index_dir}/out/other.fastq $file3
         fi
         
         mv {params.index_dir}/out/aligned.log {log} 
         rm -f {params.index_dir}/out/aligned.fastq
         rm -f {params.index_dir}/out/other.fastq
         """


def check_filtering(wildcards):
    log_file = f"{BASE}/logs/{wildcards.SRR_ID}.filtered.log"
    # Check filtered files are available
    base = f"{BASE}/filtered/{wildcards.SRR_ID}"
    if os.path.exists(base):
        files = glob_wildcards(os.path.join(base, "{name}.fastq")).name
        if len(files) == 0:
            # Folder is empty and remove it
            os.remove(log_file)
            print("Log file deleted to force filtering")
    else:
        try:
            os.remove(log_file)
        except FileNotFoundError:
            pass

    return log_file


def method_reads(wildcards):
    if config["filter_rrna"] == 1:
        return expand("{folder}/filtered/{srr}/{file}.fastq", folder=BASE,
                      srr=wildcards.SRR_ID,
                      file=glob_wildcards(os.path.join(BASE, "filtered/{"
                                                             "SRR_ID}",
                                                       "{name}.fastq")).name)
    else:
        return get_fastq(wildcards)
