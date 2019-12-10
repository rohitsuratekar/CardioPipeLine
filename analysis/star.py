#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# All functions related to STAR will go here

from constants.system import NO_OF_THREADS
from helpers.logging import Log
from helpers.resolver import PathUtil
from helpers.utils import process


def generate_index(pu: PathUtil, log: Log, force=False):
    # Make path for the index
    pu.make(pu.star_index)

    if force:
        pu.remove_all(pu.star_index)
        log.info("Old STAR index removed")
        pu.make(pu.star_index)

    opts = [
        pu.star,
        "--runMode",  # Define run mode
        "genomeGenerate",  # Mode to generate index
        "--runThreadN",  # Number of threads for running this program
        str(NO_OF_THREADS),  # If not known keep it as 1
        "--genomeDir",  # Folder where index should be made
        pu.star_index,
        "--genomeFastaFiles",  # Fasta file of full genome
        pu.genome_fasta,  # e.g. Danio_rerio.GRCz11.dna_sm.primary_assembly.fa
        "--sjdbGTFfile",  # GTF Annotation file
        pu.gtf_annotation  # e.g. Danio_rerio.GRCz11.98.chr.gtf
    ]
    log.info("STAR indexing started")
    process(opts, log,
            success=f"STAR indexing successfully performed and stored in "
                    f"{pu.star_index}",
            error="Problem occurred while generating STAR index")


def generate_alignment(sra_id: str, pu: PathUtil, log: Log):
    fastqs = pu.fastq(sra_id)
    if len(fastqs) != 2:
        log.error(f"There are total of {len(fastqs)} fastq files found in "
                  f"the default folder. Current STAR alignment need exactly "
                  f"two input fastq files. Please check your SRA id {sra_id}.")

    forward_reads = None
    reversed_reads = None
    for f in fastqs:
        if str(f).lower().endswith("1.fastq"):
            forward_reads = f
        else:
            reversed_reads = f

    if forward_reads is None or reversed_reads is None:
        log.error("Either forward or revers read has non-default name. "
                  "Forward read should end with '1.fastq' while reverse read "
                  "should end with '2.fastq'.")

    log.info("STAR mapping started on the reference genome")
    opts = [
        pu.star,
        "--runThreadN",  # Number of threads for running this program
        str(NO_OF_THREADS),  # If not known keep it as 1
        "--genomeDir",  # Path of the STAR index folder
        pu.star_index,
        "--readFilesIn",  # Input reads
        forward_reads,  # Forward read
        reversed_reads,  # Reversed read
        "--outFileNamePrefix",  # Output file prefix
        pu.star_bam_prefix(sra_id),
        "--outSAMtype",  # Direct file conversion
        "BAM",  # Directly get BAM file
        "SortedByCoordinate"  # Sort that BAM file with coordinates
    ]

    process(opts, log,
            success=f"STAR alignment is successfully completed for SRA id "
                    f"{sra_id}",
            error=f"Problem occurred while mapping {sra_id} with STAR")


def run():
    pu = PathUtil()
    log = Log(show_log=True)
    sra_id = "SRX4720626"
    generate_alignment(sra_id, pu, log)
