import os
import shutil


def get_fastq_inputs(wildcards):
    path = f"{config['srr_path']}/{wildcards.SRR_ID}/fastq"
    # Delete the fastq folder if it is empty
    if os.path.exists(path):
        if not os.listdir(path):
            shutil.rmtree(path)

    # Checkpoint should create the fastq files if there is no folder
    out_dir = checkpoints.convert_to_fastq.get(**wildcards).output.out

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


checkpoint filter_rrna:
    input:
         refs=get_rrna_reference,
         reads=get_fastq_inputs
    output:
          out=directory(config['srr_path'] + "/{SRR_ID}/filtered")

    threads: config['threads']

    shell:
         "{config[tools][sortmerna]} "
         "--split-files "  # Split files according to the reads
         "--progres "  # Show progress
         "--skip-technical "  # Skip technical reads
         "--outdir {output.out}"

rule sort_rrna:
    input: get_fastq_inputs
    output: "{SRR_ID}.sra"
    shell: "echo {wildcards.SRR_ID}"
