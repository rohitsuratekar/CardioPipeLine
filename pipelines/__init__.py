#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#

from pipelines.filter_rrna import RRNAFiltering
from pipelines.kallisto import Kallisto
from pipelines.prepare_files import PrepareRNAseq
from pipelines.salmon import Salmon
from pipelines.star import Star
from pipelines.cleanup import CleanUp
