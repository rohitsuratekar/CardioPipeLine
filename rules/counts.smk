"""
CardioPipeLine (c) 2020

Author: Rohit Suratekar, ZDG Lab, IIMCB

This deals with all rules related to generation of counts
"""

# Important: Script rules are relative to this file instead working dir.
# Hence we need to use '../scripts'

rule generate_star_counts:
    input:
         samples=config["samples"],
         gtf=config["genome"]["gtf_annotation"],
         bam="{BASE}/bams/{SRR_ID}.sra.Aligned.sortedByCoord.out.bam"
    params:
          is_paired=lambda wildcards: PIPELINE.is_paired(wildcards.SRR_ID)
    output:
          out="{BASE}/deseq2/counts/star/{SRR_ID}.star.counts"
    threads: config["threads"]
    script: "../scripts/bam_to_counts.R"

# Python Script to generate count matrix
#  {input.binary} "  #Python script 'prepDE.py'
#  -i : Input file containing gtf file locations
#  -g : Gene count Matrix
#  -t : Transcript count Matrix

rule generate_stringtie_counts:
    input:
         binary=config["tools"]["stringtie-count"],
         file="{BASE_FOLDER}/mappings/stringtie/{SRR_ID}/{SRR_ID}_assembled.gtf"
    threads: config["threads"]
    output:
          gmat="{BASE_FOLDER}/deseq2/counts/stringtie/{SRR_ID}.gene_count_matrix.csv",
          tmat="{BASE_FOLDER}/deseq2/counts/stringtie/{SRR_ID}.transcript_count_matrix.csv"
    shell:
         """
         echo "{wildcards.SRR_ID} {wildcards.BASE_FOLDER}/mappings/stringtie/{wildcards.SRR_ID}/{wildcards.SRR_ID}_assembled.gtf" > input.temp
         
         mv input.temp {wildcards.BASE_FOLDER}/mappings/stringtie/{wildcards.SRR_ID}/{wildcards.SRR_ID}_input.temp
             
         {input.binary} -i {wildcards.BASE_FOLDER}/mappings/stringtie/{wildcards.SRR_ID}/{wildcards.SRR_ID}_input.temp -g {output.gmat} -t {output.tmat}
        
         rm -f {wildcards.BASE_FOLDER}/mappings/stringtie/{wildcards.SRR_ID}/{wildcards.SRR_ID}_input.temp
         """
