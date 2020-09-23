"""
CardioPipeLine (c) 2020

Author: Rohit Suratekar, ZDG Lab, IIMCB

This deals with all rules related to Quality Control reports.
"""

# Here deliberately using temporary file so that this will run every time
rule generate_multiqc:
    input: ancient(PIPELINE.filename)
    output: temp("{BASE}/quality.html")
    shell:
         """
         multiqc {wildcards.BASE} -f -n {wildcards.BASE}/quality_report.html 
         
         touch {wildcards.BASE}/quality.html
         """
