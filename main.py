#  Copyright (c) 2020, CardioPipeLine
#
#  Author: Rohit Suratekar
#  Institute: Winata Lab, IIMCB Warsaw
#
# Main file for handling all the configuration

import subprocess

subprocess.run(["snakemake", "--core", "all"])
