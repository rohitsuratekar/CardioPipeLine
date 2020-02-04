import pandas as pd

configfile: "config/config.yaml"

wildcard_constraints:
                    SRR_ID="[SRR].*"


def get_inputs(wildcards):
    samples = pd.read_csv("config/samples.csv")
    runs = samples["Run"].values
    outputs = []
    for srr in runs:
        out = f"{str(srr).strip()}.sra"
        outputs.append(out)
    return outputs


rule run_workflow:
    input: get_inputs
    shell:
         "echo Workflow Finished"

include: "rules/preprocess.smk"
include: "rules/star.smk"
include: "rules/salmon.smk"
include: "rules/kallisto.smk"
