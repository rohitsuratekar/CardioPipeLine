import pandas as pd

configfile: 'config/config.yaml'

wildcard_constraints:
                    SRR_ID="[SRR].*"

SAMPLES_DF = pd.read_csv(config['samples'])
BASE = config['base']


def is_paired(wildcards) -> bool:
    """
    Checks if given run is with paired-end or single-end
    :return: True, if it is paired end
    """
    for _, row in SAMPLES_DF.iterrows():
        run = row["run"]
        if run.strip() == wildcards.SRR_ID:
            return row["is_paired"]
    raise KeyError(f"{wildcards.SRR_ID} not found in 'samples.csv'")


def get_final_outputs(wildcards):
    files = []
    base_f = f"{config['base']}/filtered"
    for _, row in SAMPLES_DF.iterrows():
        srr = row["run"]
        if row["is_paired"]:
            files.append(f"{base_f}/{srr}.sra.filtered_1.fastq")
            files.append(f"{base_f}/{srr}.sra.filtered_2.fastq")
        else:
            files.append(f"{base_f}/{srr}.filtered.fastq")
    return files


"""
Important: Rule order is important. Be careful with which rule you use before.
"""
rule all:
    input: get_final_outputs
    threads: config["threads"]

# Include all other rules files
include: "rules/ncbi.smk"
include: "rules/sortmerna.smk"