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

# It is important to use SRA tools 2.10.0 and above because change in Name
# resolution methods. See Issue : https://github.com/ncbi/sra-tools/issues/253
TOOL_SRA = f"{DEX_WORK}/Tools/sratoolkit.2.10.0-ubuntu64/bin"

TOOL_SORT_ME_RNA = f"{DEX_WORK}/Tools/sortmerna-4.0.0/bin"
TOOL_UNMERGE = f"{DEX_WORK}/CardioPipeLine/external/unmerge-paired-reads.sh"

TOOL_SALMON = f"{DEX_WORK}/Tools/salmon-latest_linux_x86_64/bin"

# Script needed to generate Salmon Decoy index
# https://github.com/COMBINE-lab/SalmonTools/blob/master/scripts/generateDecoyTranscriptome.sh
TOOL_SALMON_DECOY = f"{DEX_WORK}/CardioPipeLine/external/generateDecoyTranscriptome.sh"

TOOL_MASHMAP = f"{DEX_WORK}/Tools/mashmap-Linux64-v2.0"
TOOL_BEDTOOLS = f"{DEX_WORK}/Tools/bedtools2/bin"
TOOL_KALLISTO = f"{DEX_WORK}/Tools/kallisto"
