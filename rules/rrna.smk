import os
import shutil

wildcard_constraints:
                    SRR_ID="[SRR].*"

checkpoint convert_to_fastq:
    params:
          fasterq=config['tools']['fasterq-dump']
    output:
          out=directory(expand("{base}/{srr}/fastq", base=config['base'],
                               srr="{SRR_ID}"))
    threads: config['threads']
    shell:
         "{params.fasterq} {wildcards.SRR_ID} "  # Input sra ID
         "--split-files "  # Split files according to the layout 
         "--progres "  # Show progress
         "--skip-technical "  # Skip technical reads
         "--threads {threads} "  # Number of threds
         "--outdir {output.out}"
         # Output dir to store the converted fastq files


def get_fastq_inputs(wildcards):
    path = f"{config['base']}/srr/{wildcards.SRR_ID}/fastq"
    # Delete the fastq folder if it is empty
    if os.path.exists(path):
        if not os.listdir(path):
            shutil.rmtree(path)

    # Checkpoint should create the fastq files if there is no folder
    out_dir = checkpoints.convert_to_fastq.get(**wildcards).output.out

    print(f"Thois is here {out_dir}")

    return expand("{path}/{file}",
                  path=out_dir,
                  file=glob_wildcards("{name}.fastq").name)


def get_rrna_reference(wildcards):
    files = []
    db_path = config['rrna_db']
    for file in os.listdir(db_path):
        if str(file).endswith(".fasta"):
            if str(file).startswith("silva-euk"):
                files.append(f"{db_path}/{file}")
    return files


rule sort_rrna:
    input: get_fastq_inputs
    output: "{SRR_ID}.sra"
    shell: "echo {wildcards.SRR_ID} here"
