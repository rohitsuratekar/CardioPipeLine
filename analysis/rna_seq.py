#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# RNA seq analysis

from helpers import ConfigParser, MetaParser
from pipelines import *


class RNASeq:
    def __init__(self, sra_id: str, config: ConfigParser):
        self.sra = sra_id.strip()
        self.config = config
        self.meta = MetaParser(
            self.sra,
            config.names.sra.folder(self.sra),
            config.log)

        self.salmon = True
        self.star = True
        self.kallisto = True
        self.cleanup = True

    def run(self):
        self.config.log.info(f"RNA-seq analysis of {self.sra} started")
        PrepareRNAseq(self.config, self.meta).run()
        RRNAFiltering(self.config, self.meta).run()
        if self.salmon:
            Salmon(self.config, self.meta).run()
        if self.kallisto:
            Kallisto(self.config, self.meta).run()
        if self.star:
            Star(self.config, self.meta).run()
            StringTie(self.config, self.meta).run()
        if self.cleanup:
            CleanUp(self.config, self.meta).run()
        self.config.log.info(f"RNA-seq analysis of {self.sra} finished")
