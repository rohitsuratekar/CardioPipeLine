import pandas as pd

configfile: "config/config.yaml"

wildcard_constraints:
                    SRR_ID="[SRR].*"

SAMPLES = pd.read_csv("config/samples.csv")


def layout_expand(srr: str):
    data = dict(zip(SAMPLES["Run"], SAMPLES["LibraryLayout"]))
    return data[srr.strip()] == "PAIRED"


def output_file(method: str, srr: str):
    # Sanity check
    if method not in ["star", "salmon", "kallisto", "stringtie"]:
        raise Exception(f"Method {method} is not available in this workflow")
    base = f"{config['base']}/methods/{method.strip().lower()}/{srr}"
    if method == "salmon":
        return f"{base}/quant.sf"
    elif method == "kallisto":
        return f"{base}/abundance.tsv"
    elif method == "star":
        return f"{base}/{srr.strip()}_Aligned.sortedByCoord.out.bam"
    elif method == "stringtie":
        return f"{base}/{srr.strip()}_gene_expression.tsv"


def get_inputs(wildcards):
    runs = SAMPLES["Run"].values
    outputs = []
    for srr in runs:
        outputs.append(output_file("salmon", srr))
    return outputs


rule run_workflow:
    input: get_inputs
    shell:
         "echo Workflow Finished"

include: "rules/preprocess.smk"
include: "rules/star.smk"
include: "rules/salmon.smk"
include: "rules/kallisto.smk"
