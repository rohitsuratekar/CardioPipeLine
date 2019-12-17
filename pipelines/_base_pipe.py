#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# Base Pipeline

from helpers import ConfigParser, MetaParser, Log


class PipeLine:
    def __init__(self):
        self._config = None
        self._meta = None
        self._log = None

    @property
    def config(self) -> ConfigParser:
        return self._config

    @property
    def meta(self) -> MetaParser:
        return self._meta

    @property
    def log(self) -> Log:
        if self._log is None:
            self._log = Log(show_log=True)
        return self._log

    @config.setter
    def config(self, value: ConfigParser):
        self._config = value
        if self._log is None:
            self._log = value.log

    @meta.setter
    def meta(self, value: MetaParser):
        self._meta = value

    def fasta_for_analysis(self, srr):
        # Check if filtered sequence is available
        filtered = self.meta.get_filtered(srr)
        if len(filtered) > 0:
            return filtered

        # Else returned unfiltered fasta
        return self.meta.get_fastq(srr)
