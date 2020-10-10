"""
CardioPipeLine (c) 2020

Author: Rohit Suratekar, ZDG Lab, IIMCB

This deals with all rules related StringTie
"""

# Run StringTie
rule run_stringtie:
    input:
         stringtie=config["tools"]["stringtie"],
         bam="{BASE_FOLDER}/bams/{SRR_ID}.sra.Aligned.sortedByCoord.out.bam",
         gtf=ancient(config["genome"]["gtf_annotation"])
    threads: config["threads"]
    output:
          out="{BASE_FOLDER}/mappings/stringtie/{SRR_ID}/{SRR_ID}_gene_expression.tsv",
          out_gtf="{BASE_FOLDER}/mappings/stringtie/{SRR_ID}/{SRR_ID}_assembled.gtf"
    shell:
         "{input.stringtie} "  # Binary
         "-p {threads} "  # Number of threads
         "-G {input.gtf} "  # Genome gtf file
         "-e "  # only estimate the abundance of given reference transcripts
         "-B "  # enable output of Ballgown table files
         "-o {output.out_gtf} "  # Output GTF file name
         "-A {output.out} "  # Gene abundances file
         "{input.bam}"  # Input BAM file
