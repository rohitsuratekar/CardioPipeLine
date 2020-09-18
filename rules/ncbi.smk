"""
CardioPipeLine (c) 2020

Author: Rohit Suratekar, ZDG Lab, IIMCB

This deals with all rules related to NCBI-tools
"""

"""
Downloads SRA files with prefetch.

If you are having difficulty in downloading SRA files with prefetch, 
download them manually and save it in your "base/sra" folder with "run_id.sra" 
as a file name. This should be same run ID which you will provide in the 
samples.csv file

"""
rule get_sra:
    input: ancient("config/config.yaml")
    threads: config["threads"]
    output:"{BASE_FOLDER}/sra/{SSR_ID}.sra"
    params:
          prefetch=config['tools']['prefetch']
    shell:
         "mkdir -p {wildcards.BASE_FOLDER}/sra && "  # Make folder if not present
         "{params.prefetch} -p {wildcards.SSR_ID} "  # Prefetch main command
         "--output-file {output} "  # Write to specific file

rule get_fastq_paired:
    input: ancient("{BASE_FOLDER}/sra/{SRR_ID}.sra")
    threads: config["threads"]
    output:"{BASE_FOLDER}/fastq/{SRR_ID}.sra_1.fastq",
          "{BASE_FOLDER}/fastq/{SRR_ID}.sra_2.fastq"
    params:
          fasterq=config["tools"]["fasterq-dump"]

    shell:
         "{params.fasterq} "  # Binary
         "{input} "  # Input SRA File
         "--split-files "  # Split files according to the reads
         "--progress "  # Show progress
         "--skip-technical "  # Skip technical reads
         "--threads {threads} "  # Number of threads
         "--outdir {wildcards.BASE_FOLDER}/fastq"

rule get_fastq_single:
    input: ancient("{BASE_FOLDER}/sra/{SRR_ID}.sra")
    threads: config["threads"]
    output: "{BASE_FOLDER}/fastq/{SRR_ID}.sra.fastq"
    params:
          fasterq=config["tools"]["fasterq-dump"]
    shell:
         "{params.fasterq} "  # Binary
         "{input} "  # Input SRA File
         "--split-files "  # Split files according to the reads
         "--progress "  # Show progress
         "--skip-technical "  # Skip technical reads
         "--threads {threads} "  # Number of threads
         "--outdir {wildcards.BASE_FOLDER}/fastq"
