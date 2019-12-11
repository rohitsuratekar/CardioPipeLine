#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# All functions related to STAR will go here

import os

from constants.system import NO_OF_THREADS
from helpers.logging import Log
from helpers.resolver import PathUtil
from helpers.utils import process


def check_if_index_exist(pu: PathUtil, log: Log):
    if pu.exists(pu.star_index):
        files = os.listdir(pu.star_index)
        if len(files) > 0:
            log.info("STAR index folder is non empty and assuming it has "
                     "index files. If you encounter any problem while using "
                     "this index, use 'force=True' option to recreate it.")
            return True
        else:
            log.info(f"Default STAR index folder is empty {pu.star_index}")
            return False
    else:
        log.info(f"STAR index does not exist in path {pu.star_index}")
        return False


def generate_index(pu: PathUtil, log: Log, force=False):
    # Check if index files already exists

    if check_if_index_exist(pu, log) and not force:
        log.info("Skipping STAR index generation")
        return

        # Make path for the index
    pu.make(pu.star_index)

    if force:
        log.info("'Force=True' is provided and hence recreating of STAR index")
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
        pu.genome_fasta,  # e.g. Danio_rerio.GRCz11.dna.primary_assembly.fa
        "--sjdbGTFfile",  # GTF Annotation file
        pu.gtf_annotation  # e.g. Danio_rerio.GRCz11.98.chr.gtf
    ]
    log.info("STAR indexing started")
    process(opts, log,
            success=f"STAR indexing successfully performed and stored in "
                    f"{pu.star_index}",
            error="Problem occurred while generating STAR index")

    log.info("STAR generates 'Log.out' as a log file which can be found in "
             "the root of this script")


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
    generate_index(pu, log)
