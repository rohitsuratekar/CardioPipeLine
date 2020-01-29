import pandas as pd

include: "rules/common.smk"


def run_all():
    samples = pd.read_csv("config/samples.csv")
    expected_outputs = []
    for sample in samples.itertuples():
        srr = getattr(sample, "Run")
        path = config["srr_path"]
        expected_outputs.append(f"{srr}.sra")
    return expected_outputs


rule run_workflow:
    input: run_all()

include: "rules/ncbi.smk"
include: "rules/rrna.smk"
