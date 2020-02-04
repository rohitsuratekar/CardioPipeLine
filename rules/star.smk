BASE = config['base']

rule make_star_index:
    input:
         genome=config['genome']['fasta'],
         gtf=config['genome']['gtf_annotation'],
         star=config["tools"]["star"],

    output: expand("{folder}/index/star/SAindex", folder=BASE)

    params: folder="{}/index/star".format(BASE)

    threads: config['threads']

    shell:
         "{input.star} --runMode genomeGenerate "  # STAR genome mode
         "--runThreadN {threads} "  # Threads
         "--genomeDir {params.folder} "  # Folder where index will be generated
         "--genomeFastaFiles {input.genome} "  # Genome fasta file
         "--sjdbGTFfile {input.gtf}"  # GTF Annotation file

rule run_star:
    input:
         star=config["tools"]["star"],
         star_index=expand("{folder}/index/star/SAindex", folder=BASE),
         files=method_input_files
    output:
          bam=expand("{folder}/methods/star/{srr}/{"
                     "srr}_Aligned.sortedByCoord.out.bam",
                     folder=BASE,
                     srr="{SRR_ID}")

    params:
          prefix=expand("{folder}/methods/star/{srr}/{srr}_",
                        folder=BASE,
                        srr="{SRR_ID}")

    threads: config["threads"]

    shell:
         "{input.star} "  # Binary
         "--runThreadN {threads} "  # Number of threads
         "--genomeDir {BASE}/index/star "  # Path to index folder
         "--outFileNamePrefix {params.prefix} "  # Prefix of the ouput file
         "--outSAMtype BAM "  # Directly output BAM file
         "SortedByCoordinate "  # Sort BAM by coordinates
         "--readFilesIn {input.files}"  # Input reads

rule run_stringtie:
    input:
         stringtie=config["tools"]["stringtie"],
         bam=expand("{folder}/methods/star/{srr}/{"
                    "srr}_Aligned.sortedByCoord.out.bam",
                    folder=BASE,
                    srr="{SRR_ID}"),
         gtf=config["genome"]["gtf_annotation"]
    threads: config["threads"]
    output:
          out=expand("{folder}/methods/stringtie/{srr}/{"
                     "srr}_gene_expression.tsv",
                     folder=BASE, srr="{SRR_ID}"),
          out_gtf=expand("{folder}/methods/stringtie/{srr}/{"
                         "srr}_assembled.gtf",
                         folder=BASE, srr="{SRR_ID}")
    shell:
         "{input.stringtie} "  # Binary
         "-p {threads} "  # Number of threads
         "-G {input.gtf} "  # Genome gtf file
         "-e "  # only estimate the abundance of given reference transcripts
         "-B "  # enable output of Ballgown table files
         "-o {output.out_gtf} "  # Output GTF file name
         "-A {output.out} "  # Gene abundances file
         "{input.bam}"  # Input BAM file
