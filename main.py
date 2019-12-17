#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# Main file to control the project

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

    all_ids = ["SRX4720625", "SRX4720626",
               "SRX4720628", "SRX4720629"]

    for sra in all_ids:
        RNASeq(sra, config).run()
