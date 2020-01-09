#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# Differential analysis with DESeq2

from typing import List

from helpers import ConfigParser, MetaParser
from pipelines import CountMatrix, DESeqMeta


class DESeq2:
    def __init__(self, config: ConfigParser, method: str):
        self.config = config
        self.samples = []
        self._conditions = {}

        self.method = method.strip().lower()

        self.salmon = False
        self.star = False
        self.kallisto = False
        self.stringtie = False

        self.log = config.log

        self._assign_method()

    def add_samples(self, sample_list: List[str]):
        self.samples = sample_list

    def add_condition(self, name: str, condition_list: List[str]):
        self._conditions[name] = condition_list

    def _assign_method(self):
        if self.method == "star":
            self.star = True
        elif self.method == "salmon":
            self.salmon = True
        elif self.method == "kallisto":
            self.kallisto = True
        elif self.method == "stringtie":
            self.stringtie = True
        else:
            self.config.log.error(f"Method {self.method} is not supported in "
                                  f"the DESeq2 pipeline yet. Available "
                                  f"options are [star, salmon, kallisto, "
                                  f"stringtie]")

    def _generate_count_matrix(self, sra):

        meta = MetaParser(
            sra,
            self.config.names.sra.folder(sra),
            self.config.log)

        if self.salmon:
            CountMatrix(self.config, meta, "salmon").run()
        if self.kallisto:
            CountMatrix(self.config, meta, "kallisto").run()
        if self.star:
            CountMatrix(self.config, meta, "star").run()
        if self.stringtie:
            CountMatrix(self.config, meta, "stringtie").run()

    def _generate_meta_files(self):
        deq = DESeqMeta(self.config, self.method)
        deq.samples(self.samples)

        for condition in self._conditions:
            deq.add_conditions(condition, self._conditions[condition])

        deq.run()

    def _check_setup(self):
        samples = len(self.samples)
        if samples == 0:
            self.log.error("No sample is provided for the analysis. Please "
                           "add samples before running the pipeline")

        if len(self._conditions) == 0:
            self.log.error("No conditions are provided for the analysis. "
                           "Please add at least one condition to run the "
                           "pipeline")

        for k in self._conditions.values():
            if len(k) != samples:
                self.log.error("Number of 'condition-details' should be same "
                               "as number of samples.")

    def run(self):
        self.config.log.info(f"DESeq2 analysis started")

        # Check setup
        self._check_setup()

        # Generate individual count matrix
        for sra in self.samples:
            self._generate_count_matrix(sra)

        # Generate merged count matrix and meta files
        self._generate_meta_files()

        self.config.log.info(f"DESeq2 analysis finished")
