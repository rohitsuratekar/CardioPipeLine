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


def sort_me_rna(sra_id: str, pu: PathUtil, log: Log):
    opts = [
        pu.sortmerna  # SortMeRNA Binary path
    ]
    # Add references
    for f in os.listdir(pu.rrna_db):
        if str(f).endswith(".fasta"):
            opts.append("-ref")  # Option to add rrna reference file
            opts.append(f"{pu.rrna_db}/{f}")  # Here you can filter rrna
            # files if you want

    # Add input files
    file_count = 0
    for f in pu.fastq(sra_id):
        opts.append("-reads")  # Option to add reads
        opts.append(f)  # Upstream Fastq files
        file_count += 1

    # Sanity check
    if file_count not in [1, 2]:
        log.error(f"rRNA sorting is not supported for current sample as it "
                  f"has total of {file_count} fastq files. Please provide 1 "
                  f"or 2 file for single and paired end sequencing.")

    # Add other options
    opts.extend([
        "-workdir",  # Working folder for storing index and output
        pu.sortmerna_index,
        "-fastx",  # Output reads into fastq format
        "-a",  # Number of threads
        str(NO_OF_THREADS),  # If not available, use 1
        "-v",  # Produce verbose output
        "-paired_in"  # This will put all non-rrna reads in other.fasta
    ])

    process(opts, log,
            success="rRNA filtering with SortMeRNA is successful",
            error="Something went wrong in rRNA filtering with SortMeRNA")


def run():
    pu = PathUtil()
    log = Log(show_log=True)
    sra_id = "SRX4720626"
    sort_me_rna(sra_id, pu, log)
