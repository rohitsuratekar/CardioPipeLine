#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# All the constants related to tools
#
# Two constants DEX_WORK and DEX_DATA are taken from local.py files which
# are private. You should change these paths here according to your system.
#
# Avoid trailing slash in the path. It should be handled at the usage level

from local import *

TOOL_STAR = f"{DEX_WORK}/Tools/STAR-2.7.3a/bin/Linux_x86_64"
TOOL_SRA = f"{DEX_WORK}/Tools/sratoolkit.2.10.0-ubuntu64/bin"
