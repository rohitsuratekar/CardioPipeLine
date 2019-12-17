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

    sra = "SRX4720625"

    RNASeq(sra, config).run()
