#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# rRNA filtering related functions

import os

from constants.system import NO_OF_THREADS
from helpers.logging import Log
from helpers.resolver import PathUtil
from helpers.utils import process


def post_process(sra_id: str, pu: PathUtil, log: Log):
    # Sanity check
    if not pu.exists(pu.files.sortmerna_out(pu.index.sortmerna)):
        log.error(f"SortMeRNA postprocessing stopped because output files are "
                  f"not found in the default folder {pu.index.sortmerna}")

    opts = [
        pu.tools.unmerge,  # Path to un-merge script
        pu.files.sortmerna_out(pu.index.sortmerna),  # Merge file
        pu.files.filtered_read(1, sra_id, True),  # Read 1
        pu.files.filtered_read(2, sra_id, True)  # Read 2
    ]

    process(opts, log, success="Merged files successfully unmerged",
            error="Something went wrong in un-merging files")

    # Move log file to SRA folder
    log.info("Moving log file to filtered folder")
    pu.move(pu.files.sortmerna_log(pu.index.sortmerna),
            pu.files.filtering_log(sra_id))

    # Delete filtered files
    pu.remove(pu.files.sortmerna_out(pu.index.sortmerna))
    pu.remove(pu.files.sortmerna_align(pu.index.sortmerna))

    log.info("Temporary filtered files removed")


def sort_me_rna(sra_id: str, pu: PathUtil, log: Log):
    # Check if indexing is already done

    opts = [
        pu.tools.sortmerna  # SortMeRNA Binary path
    ]
    # Add references
    for f in os.listdir(pu.rrna_db):
        if str(f).endswith(".fasta"):
            # For faster filtering, we will only use Eukaryotic rRNA data
            if str(f).startswith("silva-euk"):
                opts.append("-ref")  # Option to add rrna reference file
                opts.append(f"{pu.rrna_db}/{f}")  # Here you can filter rrna
                # files if you want

    # Add input files
    fastq_files = sorted(
        pu.fastq(sra_id))  # Sort it so that forward read goes first
    for f in fastq_files:
        opts.append("-reads")  # Option to add reads
        opts.append(f)  # Upstream Fastq files

    # Sanity check
    if len(fastq_files) not in [1, 2]:
        log.error(f"rRNA sorting is not supported for current sample as it "
                  f"has total of {len(fastq_files)} fastq files. Please"
                  f" provide 1 or 2 file for single and paired end "
                  f"sequencing.")

    # Add other options
    # if "-paired_in" is used, both reads will be merged, hence don't use it
    opts.extend([
        "-workdir",  # Working folder for storing index and output
        pu.index.sortmerna,
        "-fastx",  # Output reads into fastq format
        "-other",  # Output the Other.fastq file (which has NON-rRNA reads)
        "-a",  # Number of threads
        str(NO_OF_THREADS),  # If not available, use 1
        "-v",  # Produce verbose output
        "-num_alignments",  # Only first alignment passing of E value
        str(1),  # This will be fastest choice when only filtering is needed,
        "-paired_in"
        # if either of read match, put both reads into the 'aligned.fasta'
        # file.
    ])

    # Delete the kvdb folder if exists. This is needed for SortMeRNA
    pu.remove(f"{pu.index.sortmerna}/kvdb")

    process(opts, log,
            success="rRNA filtering with SortMeRNA is successful",
            error="Something went wrong in rRNA filtering with SortMeRNA")


def run():
    pu = PathUtil()
    log = Log(show_log=True)
    sra_id = "SRX4720626"
    post_process(sra_id, pu, log)
