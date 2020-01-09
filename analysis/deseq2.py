#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# Differential analysis with DESeq2

import subprocess
from typing import List

from helpers import ConfigParser, MetaParser
from pipelines import CountMatrix, DESeqMeta


class DESeq2:
    def __init__(self, config: ConfigParser, method: str):
        self.config = config
        self.samples = []
        self._conditions = {}
        self.design = None
        self.deseq_test = "Wald"
        self.reduced_design = "~ 1"
        self.contrast = None
        self.method = method.strip().lower()
        self.alpha = 0.1

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

    def add_design(self, design: str):
        if len(design.strip()) == 0:
            self.log.error("Design can not be empty string")
        self.design = design

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

        if self.design is None:
            self.log.error("No design formula is provided. Please add design "
                           "formula for the DESeq2 analysis.")

    def _perform_analysis(self):
        opts = [
            "Rscript",
            "rscripts/deseq2.R",
            "-c",
            f"{self.config.deseq2_final_output}/{self.method}.counts.csv",
            "-m",
            f"{self.config.deseq2_final_output}/{self.method}.meta.csv",
            "-d",
            self.design,
            "-t",
            self.deseq_test,
            "-r",
            self.reduced_design,
            "-o",
            f"{self.config.deseq2_final_output}/{self.method}"
        ]

        if self.contrast is not None:
            if len(self.contrast) != 3:
                self.log.error("Contrast should contain exactly 3 string "
                               "elements")
            opts.extend([
                "-n",
                "\t".join(self.contrast),
                "-a",
                str(self.alpha)
            ])

        if subprocess.run(opts).returncode != 0:
            self.log.error("Something went wrong in performing DESeq2 "
                           "analysis")

        self.log.info(f"Successful DESeq2 analysis. All the results are "
                      f"stored in {self.config.deseq2_final_output}")

    def run(self):
        self.config.log.info(f"DESeq2 analysis started")

        # Check setup
        self._check_setup()

        # Generate individual count matrix
        for sra in self.samples:
            self._generate_count_matrix(sra)

        # Generate merged count matrix and meta files
        self._generate_meta_files()

        # Perform the analysis
        self._perform_analysis()

        self.config.log.info(f"DESeq2 analysis finished")
