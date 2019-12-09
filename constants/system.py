#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# All system level constants
#
# Two constants DEX_HOME, DEX_WORK and DEX_DATA are taken from local.py files
# which are private. You should change these paths here according to your
# system.
#
# Avoid trailing slash in the path. It should be handled at the usage level

from local import *

SHOW_LOG = True
SUCCESS = 0

FOLDER_NCBI_PUBLIC = f"{DEX_HOME}/ncbi/public"
FOLDER_SRA = f"{DEX_DATA}/sra"
