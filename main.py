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

from analysis.rna_seq import RNASeq
from helpers import ConfigParser

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

    for sra in all_ids:
        RNASeq(sra, config).run()
