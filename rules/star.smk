"""
CardioPipeLine (c) 2020

Author: Rohit Suratekar, ZDG Lab, IIMCB

This deals with all rules related STAR analysis
"""

BASE = config['base']

# Generate the STAR index
rule make_star_index:
    input:
         genome=ancient(config['genome']['fasta']),
         gtf=ancient(config['genome']['gtf_annotation']),
         star=ancient(config["tools"]["star"]),

    output: expand("{folder}/index/star/SAindex", folder=BASE)

    params: folder="{}/index/star".format(BASE)

    threads: config['threads']

    shell:
         "{input.star} --runMode genomeGenerate "  # STAR genome mode
         "--runThreadN {threads} "  # Threads
         "--genomeDir {params.folder} "  # Folder where index will be generated
         "--genomeFastaFiles {input.genome} "  # Genome fasta file
         "--sjdbGTFfile {input.gtf}"  # GTF Annotation file

# Run the alignment
# rule run_star:
#     input:
#          star=config["tools"]["star"],
#          star_index=ancient("{base_folder}/index/star/SAindex"),
#          files=method_reads,
#     output:
#           bam="{base_folder}/bams/{SRR_ID}.sra.Aligned.sortedByCoord.out.bam"
#
#     params:
#           prefix=expand("{folder}/bams/{srr}.sra.", folder=BASE,
#                         srr="{SRR_ID}")
#
#     threads: config["threads"]
#
#     shell:
#          "{input.star} "  # Binary
#          "--runThreadN {threads} "  # Number of threads
#          "--genomeDir {BASE}/index/star "  # Path to index folder
#          "--outFileNamePrefix {params.prefix} "  # Prefix of the ouput file
#          "--outSAMtype BAM "  # Directly output BAM file
#          "SortedByCoordinate "  # Sort BAM by coordinates
#          "--readFilesIn {input.files}"  # Input reads
