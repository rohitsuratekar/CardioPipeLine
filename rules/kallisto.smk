BASE = config['base']

rule index_kallisto:
    input:
         kallisto=config["tools"]["kallisto"],
         transcript=config["genome"]["transcript"]

    output: f"{BASE}/index/kallisto/kallisto.idx"

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
        return f"--single -l {frg_length} -s {frg_sd} {files[0]}"


rule kallisto_quant:
    input:
         kallisto=config["tools"]["kallisto"],
         k_index=f"{BASE}/index/kallisto/kallisto.idx",
         gtf=config["genome"]["gtf_annotation"],
         files=method_input_files
    threads: config["threads"]
    output:
          expand("{folder}/methods/kallisto/{srr}/abundance.tsv",
                 folder=BASE, srr="{SRR_ID}")
    params:
          reads=lambda wildcards, input: make_kallisto_reads(input.files)
    shell:
         "{input.kallisto} quant "  # Kallisto quant mode
         "-i {input.k_index} "  # Path to the index file
         "-o {BASE}/methods/kallisto/{wildcards.SRR_ID} "  # Output folder
         "-g {input.gtf} "  # Genome GTF file
         "-t {threads} "  # Number of threads
         "--bias "  # Sequence based bias correction,
         "--plaintext "  # Plaintext output instead HDF5
         "{params.reads}"  # Add reads and layout specific options

rule dry_run:
    input:
         expand("{folder}/methods/kallisto/{srr}/abundance.tsv",
                folder=BASE, srr="{SRR_ID}")
    output: "{SRR_ID}.sra"
    shell:
         "touch {wildcards.SRR_ID}.sra"
