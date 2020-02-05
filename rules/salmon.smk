BASE = config['base']

rule generate_decoy:
    input:
         gtf=config["genome"]["gtf_annotation"],
         fasta=config["genome"]["fasta"],
         transcript=config["genome"]["transcript"],
         bedtools=config["tools"]["bedtools"],
         mashmap=config["tools"]["mashmap"]
    output:
          out_index=f"{BASE}/index/salmon/gentrome.fa",
          out_keys=f"{BASE}/index/salmon/decoys.txt"
    threads: config["threads"]
    shell:
         "scripts/generateDecoyTranscriptome.sh "
         "-j {threads} "  # Number of threads
         "-a {input.gtf} "  # Genome GTF File
         "-g {input.fasta} "  # Genome in fasta
         "-t {input.transcript} "  # Transcript file in fasta
         "-o {BASE}/index/salmon "  # Salmon index output folder
         "-b {input.bedtools} "  # Path to bedtools
         "-m {input.mashmap} "  # Path to mashmap

rule index_salmon:
    input:
         decoy_index=f"{BASE}/index/salmon/gentrome.fa",
         decoy_keys=f"{BASE}/index/salmon/decoys.txt",
         salmon=config["tools"]["salmon"]
    threads: config["threads"]
    output: f"{BASE}/index/salmon/pos.bin"
    shell:
         "{input.salmon} index "  # Salmon index mode
         "-t {input.decoy_index} "  # Decoy transcript
         "-i {BASE}/index/salmon "  # Path to index folder
         "--decoys {input.decoy_keys} "  # Decoy keys
         "-k 31 "  # K-mer (31) is recommended by the authors
         "-p {threads}"  # Number of threads

rule salmon_quant:
    input:
         salmon=config["tools"]["salmon"],
         salmon_index=f"{BASE}/index/salmon/pos.bin",  # This is just needed
         # to ensure that index is prepared.
         gtf=config["genome"]["gtf_annotation"],
         files=method_reads,
    output:
          expand("{folder}/methods/salmon/{srr}/quant.sf",
                 folder=BASE, srr="{SRR_ID}")
    params:
          reads=lambda wildcards, input: " ".join([
              "-{} {}".format(i + 1, x)
              if len(input.files) == 2 else "-r {}".format(x)
              for i, x in enumerate(sorted(input.files))
          ])
    threads: config["threads"]
    shell:
         "{input.salmon} quant "  # Salmon quant mode
         "-i {BASE}/index/salmon "  # Path of the index
         "-l A "  # Automatic library type
         "--validateMapping "  #enables selective alignment of the 
         # sequencing reads
         "--gcBias "  # Corrects for  fragment-level GC biases. Also 
         # recommended by the DESEq2 analysis
         "-p {threads} "  # Number of threads   
         "-g {input.gtf} "  # Additional GTF file to get quants.genes.sf
         "-o {BASE}/methods/salmon/{wildcards.SRR_ID} "  # output folder
         "{params.reads}"  # Input reads
