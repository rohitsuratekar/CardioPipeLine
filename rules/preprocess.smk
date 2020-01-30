checkpoint fasterq_download:
    output: directory(config['base'] + "/srr/{SRR_ID}/fastq")
    shell:
         "{config[tools][fasterq-dump]} {wildcards.SRR_ID} "
         "--split-files "  # Split files according to the reads
         "--progres "  # Show progress
         "--skip-technical "  # Skip technical reads
         "--outdir {output}"


def collect_fastq(wildcards):
    path = checkpoints.fasterq_download.get(**wildcards).output[0]
    return expand("{path}/{file}",
                  path=path,
                  file=glob_wildcards("{name}.fastq").name)


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

checkpoint filter_rrna:
    input:
         reads=collect_fastq,
         refs=expand("{path}/{files}",
                     path=config['rrna_db'],
                     files=glob_wildcards("{name}.fasta").name)
    output: directory(config['base'] + "/srr/{SRR_ID}/filtered")
    threads: config['threads']
    params:
          ref_str=lambda _, input: " ".join(
              ["-ref {}".format(x) for x in input.refs]),
          read_str=lambda _, input: " ".join(
              ["-reads {}".format(x) for x in input.reads]),
          index_dir="{}/index/sortmerna".format(config['base'])
    shell:
         """
         rm -rf {params.index_dir}/kvdb
         
         {config[tool][sortmerna]} {params.ref_str} {params.read_str} \
         -workdir {params.index_dir} -fastx -other -a {threads} \
         -num_alignment 1 -paired_in
         """

rule try_run:
    input: collect_fastq
    output: "{SRR_ID}.sra"
    shell:
         "touch {wildcards.SRR_ID}.sra"
