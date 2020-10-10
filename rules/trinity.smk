"""
CardioPipeLine (c) 2020

Author: Rohit Suratekar, ZDG Lab, IIMCB

This deals with all rules related Trinity
"""

"""
Trinity Options
--seqType : Sequence Type
--left: Left reads
--right: Right reads
--CPU: therads
--output : Output folder
--max_memory : Maximum Memory
"""

rule de_novo_trinity_paired:
    input:
         files=get_fastq_files,
         trinity=config["tools"]["trinity"],
         jellyfish=config["tools"]["jellyfish"],
         samtools=config["tools"]["samtools"],
         bowtie2=config["tools"]["bowtie2"],
         salmon=config["tools"]["salmon"]
    output: "{BASE_FOLDER}/mappings/trinity/trinity_{SRR_ID}/Trinity.fasta"
    threads: config["threads"]
    shell:
         """ 
         export PATH="$(dirname {input.samtools}):$PATH"
         export PATH="$(dirname {input.jellyfish}):$PATH"
         export PATH="$(dirname {input.bowtie2}):$PATH"
         export PATH="$(dirname {input.salmon}):$PATH"
         
         {input.trinity} --seqType fq --left {input.files[0]} --right {input.files[1]} --CPU {threads} --output {wildcards.BASE_FOLDER}/mappings/trinity/trinity_{wildcards.SRR_ID} --max_memory 50G --verbose
         """

"""
Build BowTie indexs

bowtie2-build assembled transcript output_prefix
"""
rule bowtie2_index:
    input:
         transcripts="{BASE_FOLDER}/mappings/trinity/trinity_{SRR_ID}/Trinity.fasta",
         bowtie2=config["tools"]["bowtie2-build"]
    threads: config["threads"]
    output:
          "{BASE_FOLDER}/mappings/trinity/trinity_{SRR_ID}/Trinity.fasta.1.bt2"
    shell:
         """
         {input.bowtie2} {input.transcripts} {input.transcripts}
         """

"""
Get statistics from paired-end
Bowtie2:
-p : number of threads
-q : input files are in fastq format
--no-unal: suppress SAM records for unaligned reads
-x: Index files prefix (without X.bt2 extension)

"""
rule bowtie2_pair_end:
    input:
         files=get_fastq_files,
         bowtie=config["tools"]["bowtie2"],
         samtools=config["tools"]["samtools"],
         transcripts="{BASE_FOLDER}/mappings/trinity/trinity_{SRR_ID}/Trinity.fasta",
         index="{BASE_FOLDER}/mappings/trinity/trinity_{SRR_ID}/Trinity.fasta.1.bt2"
    threads: config["threads"]
    output: "{BASE_FOLDER}/mappings/trinity/trinity_{SRR_ID}/align_stats.txt"
    shell:
         """
         {input.bowtie} -p {threads} -q --no-unal -x {input.transcripts} -1 {input.files[0]} -2 {input.files[1]} > {output}
         """
