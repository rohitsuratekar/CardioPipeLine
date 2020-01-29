wildcard_constraints:
                    SRR_ID="[SRR].*"

rule prefetch_download:
    params:
          folder=f"{config['srr_path']}"
    output:
          "{folder}/{SRR_ID}/{SRR_ID}.sra"
    shell:
         '{config[tools][prefetch]} {wildcards.SRR_ID} '
         '-O {wildcards.folder}/{wildcards.SRR_ID}'  # Define output folder

checkpoint convert_to_fastq:
    input:
         expand("{path}/{srr}/{srr}.sra",
                path=config['srr_path'],
                srr="{SRR_ID}")
    output:
          out=directory(config['srr_path'] + "/{SRR_ID}/fastq")

    threads: config['threads']

    shell:
         "{config[tools][fasterq-dump]} {input} "  #Input files
         "--split-files "  # Split files according to the reads
         "--progres "  # Show progress
         "--skip-technical "  # Skip technical reads
         "--outdir {output.out}"
