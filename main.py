#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# Main file to control the project
#
# Important Note: Make sure that external/generateDecoyTranscriptome.sh and
# external/unmerge-paired-reads.sh has execute access.

import itertools
import json
import subprocess
import sys
from collections import defaultdict

from analysis import RNASeq, DESeq2
from helpers import ConfigParser, make_path, move_path, delete_path


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


def public_deseq2_example(mtd: str,
                          sample_list: list,
                          condition_details: list,
                          contrast: list):
    config.log.info(f"Starting analysis with following parameters")
    config.log.info(f"Samples : {sample_list}")
    config.log.info(f"Conditions : {condition_details}")
    config.log.info(f"Contrast : {contrast}")

    # List of SRA ids for DESeq2 analysis
    # E.g.
    # sample_list = ["SRX3187828", "SRX3187834", "SRX3187822", "SRX3187825",
    #                "SRX3187837", "SRX3187839"]

    # Their conditions
    # In this example, I am using samples taken at different time
    # developmental time points as a condition
    # Here we are using 24hpf as a control or starting point. And by
    # convention in DESeq2 program, control is provided at the end.
    condition = "time"

    # E.g condition_details = ["48hpf", "48hpf", "48hpf", "30hpf", "30hpf",
    # "30hpf"]

    # Generate DESeq2 object with desired method
    # Currently 4 methods are available [star, salmon, kallisto, stringtie]
    deq = DESeq2(config, mtd)

    # Add samples
    deq.add_samples(sample_list)

    # Add conditions
    # Give name which you can use in statistical modelling of DESeq2 tool
    deq.add_condition(condition, condition_details)
    # You can add extra conditions like cell type, genotype extra if you
    # like. These all conditions will be available in the meta file

    # We will just do simple differential analysis based on developmental time
    design = f"~ {condition}"  # e.g "~ time"
    deq.add_design(design)

    # We can provide contrast matrix for LFC shrinkage.
    # If not given, LFC Shrinkage will not be performed
    # E.g. : ["time", "48hpf", "30hpf"]
    deq.contrast = contrast
    deq.alpha = 0.05  # Alpha usually decided p-value threshold of FDR. This
    # will affect the adjusted-p-value in the final result

    # Run the analysis
    deq.run()


def package_analysis(name: str):
    input_folder = config.deseq2_final_output
    base_path = input_folder.replace("Analysis", "DEOut")
    make_path(base_path)
    working_path = f"{base_path}/{name}"
    move_path(input_folder, working_path)
    # While using tar from another dictionary, remember to use -C option
    opts = [
        "tar",
        "-czf",
        f"{base_path}/{name}.tar.gz",
        "-C",
        base_path,
        name
    ]
    if subprocess.run(opts).returncode != 0:
        config.log.error(
            f"Something went wrong while creating archive of {working_path}")

    delete_path(working_path)
    config.log.info(f"Archive created at {working_path}")


def sample_combination():
    data = defaultdict(list)
    with open("samples.csv") as file:
        for line2 in file:
            d = line2.strip().split(',')
            data[d[1]].append(d[0])

    pairs = itertools.combinations(data.keys(), 2)
    for pair in pairs:
        current_samples = [x for x in data[pair[0]]]
        current_times = [pair[0] for _ in range(len(data[pair[0]]))]
        current_samples.extend(data[pair[1]])
        current_times.extend(
            [pair[1] for _ in range(len(data[pair[1]]))]
        )
        current_times = [f"{x}hpf" for x in current_times]
        current_contrast = sorted(pair, reverse=True)
        current_contrast = [f"{x}hpf" for x in current_contrast]
        current_contrast.insert(0, "time")
        # Make folder for the analysis
        make_path(config.deseq2_final_output)
        #
        # DESeq2 analysis
        # DESeq2 example with all 4 methods
        for method in ["star", "salmon", "kallisto", "stringtie"]:
            public_deseq2_example(method,
                                  current_samples,
                                  current_times,
                                  current_contrast)
            
        package_analysis(f"Analysis{pair[0]}v{pair[1]}")


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

    sample_combination()
