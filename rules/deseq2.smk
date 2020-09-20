"""
CardioPipeLine (c) 2020

Author: Rohit Suratekar, ZDG Lab, IIMCB

This deals with all rules related to DESeq2 analysis
"""
import pandas as pd


def get_current_input(wildcards):
    f, _, _ = PIPELINE.get_deseq2_inputs(wildcards.METHOD)
    return f


def get_current_samples(wildcards):
    _, s, _ = PIPELINE.get_deseq2_inputs(wildcards.METHOD)
    return s


def get_current_conditions(wildcards):
    _, _, c = PIPELINE.get_deseq2_inputs(wildcards.METHOD)
    return c


def generate_deseq_input(wildcards):
    files = get_current_input(wildcards)
    dfs = [pd.read_csv(x) for x in files]
    dfs = [x.set_index("gene_id") for x in dfs]
    df = pd.concat(dfs, axis=1, join="inner")
    df = df.reset_index()
    outfile = f"{config['base']}/deseq2/analysis/{wildcards.METHOD}/temp.counts"
    df.to_csv(outfile, index=False)
    return outfile


rule deseq2_with_star:
    input:
         files=get_current_input,
         combined=generate_deseq_input
    params:
          samples=lambda wildcards: get_current_samples(wildcards),
          conditions=lambda wildcards: get_current_conditions(wildcards),
          column=config['deseq2']['design_column'],
          design=config['deseq2']['design'],
          reference=config['deseq2']['reference']
    output:
          log="{BASE}/deseq2/analysis/{METHOD}/combined.counts"
    threads: config["threads"]
    script: "../scripts/deseq2.R"
