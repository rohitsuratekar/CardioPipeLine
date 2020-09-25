"""
CardioPipeLine (c) 2020

Author: Rohit Suratekar, ZDG Lab, IIMCB

This deals with all rules related STAR analysis
"""

BASE = config['base']

# Generate the STAR index
#  --runMode genomeGenerate  : STAR genome mode
#  --runThreadN : Threads
#  --genomeDir  : Folder where index will be generated
#  --genomeFastaFiles  : Genome fasta file
#  --sjdbGTFfile : GTF Annotation file
rule make_star_index:
    input:
         genome=ancient(config['genome']['fasta']),
         gtf=ancient(config['genome']['gtf_annotation']),
         star=ancient(config["tools"]["star"]),

    output: expand("{folder}/index/star/SAindex", folder=BASE)

    params: folder="{}/index/star".format(BASE)

    threads: config['threads']

    shell:
         """
         {input.star} --runMode genomeGenerate --runThreadN {threads} --genomeDir {params.folder} --genomeFastaFiles {input.genome} --sjdbGTFfile {input.gtf}
        
        if [ -f Log.out ]; then
            mkdir -p {BASE}/logs
            mv Log.out {BASE}/logs/STAR.index.log
        fi 
         """

# Run the alignment
#  --runThreadN : Number of threads
#  --genomeDir : Path to index folder
#  --outFileNamePrefix : Prefix of the ouput file
#  --outSAMtype BAM  : Directly output BAM file
#  SortedByCoordinate : Sort BAM by coordinates
#  --readFilesIn : Input reads
rule run_star:
    input:
         star=config["tools"]["star"],
         star_index=ancient("{BASE_FOLDER}/index/star/SAindex"),
         files=get_fastq_files
    output:
          bam="{BASE_FOLDER}/bams/{SRR_ID}.sra.Aligned.sortedByCoord.out.bam"

    params:
          prefix=expand("{folder}/bams/{srr}.sra.", folder=BASE,
                        srr="{SRR_ID}")

    threads: config["threads"]

    shell:
         """
         {input.star} --runThreadN {threads} --genomeDir {BASE}/index/star --outFileNamePrefix {params.prefix} --outSAMtype BAM SortedByCoordinate --readFilesIn {input.files}
         
         rm -f {BASE}/bams/{params.prefix}.Log.progress.out
         
         """
