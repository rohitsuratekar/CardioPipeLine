import pandas as pd
from models.basic import PipeLine

CONFIG_FILE = "config/config.yaml"

configfile: CONFIG_FILE

wildcard_constraints:
                    SRR_ID="[SRR].*"

SAMPLES_DF = pd.read_csv(config['samples'])
BASE = config['base']
PIPELINE = PipeLine(CONFIG_FILE)


def get_final_outputs(wildcards):
    files = []
    for _, row in SAMPLES_DF.iterrows():
        srr = row["run"]
        paired = row["is_paired"]
        files.extend(PIPELINE.output_files(srr, paired))
    return files


def get_fastq_files(wildcards):
    return PIPELINE.fastq(wildcards.SRR_ID)


"""
Important: Rule order is important. Be careful with which rule you use before.
"""
rule all:
    input: get_final_outputs
    threads: config["threads"]

# Include all other rules files
# This order is also important as some of the functions are reused in other
# files.
include: "rules/ncbi.smk"
include: "rules/sortmerna.smk"
include: "rules/star.smk"
include: "rules/stringtie.smk"
include: "rules/salmon.smk"
include: "rules/kallisto.smk"
include: "rules/counts.smk"
