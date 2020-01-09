#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# Main file to control the project
#
# Important Note: Make sure that external/generateDecoyTranscriptome.sh and
# external/unmerge-paired-reads.sh has execute access.

import json
import sys

from analysis import RNASeq, DESeq2
from helpers import ConfigParser


def public_rnaseq_example(ids: list):
    """
    This example shows how to perform RNA-seq analysis on public datasets (
    datasets which are available on NCBI's SRA database)
    https://www.ncbi.nlm.nih.gov/sra/docs/sradownload/

    You can either pass SRA id through terminal like following,
    >>>python3 main.py config.json SRX4720625
    Limitation of this method (for now) is that you can only pass single ID

    To pass multiple SRA ids, you can make 'sra.txt' file in the same folder
     and then run the program
    >>>python3 main.py config.json


    :param ids: List of SRA ids on which you want to perform RNA-seq
    analysis. e.g: SRX4720625
    """
    for sra in ids:
        # To avoid any spaces to go in the analysis and cause error
        if len(sra.strip()) > 0:
            seq = RNASeq(sra, config)
            # By default, all analysis will be performed.
            # You can selectively switch off specific analysis
            seq.star = False  # If you do not want STAR analysis

            # By default, at the end, unnecessary converted files will be
            # deleted. However, you can block this behaviour like following
            seq.cleanup = False

            # Finally run the analysis
            seq.run()


def public_deseq2_example(method: str):
    # List of SRA ids for DESeq2 analysis
    sample_list = ["SRX4720625", "SRX4720626"]

    # Their conditions
    # In this example, I am using samples taken at different time
    # developmental time points as a condition
    conditions = ["24hpf", "24hpf"]

    # Generate DESeq2 object with desired method
    # Currently 4 methods are available [star, salmon, kallisto, stringtie]
    deq = DESeq2(config, method)

    # Add samples
    deq.add_samples(sample_list)

    # Add conditions
    # Give name which you can use in statistical modelling of DESeq2 tool
    deq.add_condition("time", conditions)

    # You can add extra conditions like cell type, genotype extra if you
    # like. These all conditions will be available in the meta file

    # Run the analysis
    deq.run()


if __name__ == "__main__":
    if len(sys.argv) == 1:
        raise FileNotFoundError("config.json not found. Please pass "
                                "config.json (with full path) as an argument")

    with open(sys.argv[1]) as f:
        config = ConfigParser(json.load(f))

    all_ids = []
    if len(sys.argv) == 2:
        config.log.info("No SRA id is provided. Using 'sra.txt' file")
        with open("sra.txt") as f:
            for line in f:
                all_ids.append(str(line).strip().upper())
    else:
        all_ids.append(str(sys.argv[2]).strip().upper())
        config.log.info(f"Analysis will start for {all_ids[0]}")

    # DESeq2 example
    public_deseq2_example("star")
