configfile: "config/config.yaml"


def get_inputs(wildcards):
    srr = "SRR7882018"
    base_path = config['base']
    return [f"{srr}.sra"]


rule run_workflow:
    input: get_inputs
    shell:
         "echo Workflow Finished"

include: "rules/preprocess.smk"
