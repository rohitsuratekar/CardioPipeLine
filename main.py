#  Copyright (c) 2020, CardioPipeLine
#
#  Author: Rohit Suratekar
#  Institute: Winata Lab, IIMCB Warsaw
#
# Main file for handling all the configuration

import subprocess

subprocess.run(["snakemake", "--core", "all"])

# from models.basic import PipeLine
# import os, sys
#
# folder = os.path.dirname(sys.argv[0])
# p = PipeLine(f"{folder}/config/config.yaml")
# p.run()
