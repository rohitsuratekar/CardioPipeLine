#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# Differential analysis with DESeq2

import json
import subprocess
from typing import List

import pandas as pd

from helpers import ConfigParser, MetaParser, exists_path, make_path


class DESeq2:
    def __init__(self, config: ConfigParser, method: str):
        self.config = config
        self.samples = []
        self._conditions = {}
        self.design = None
        self.contrast = None
        self.method = method.strip().lower()
        self.alpha = 0.1

        self.log = config.log

    def add_samples(self, sample_list: List[str]):
        self.samples = sample_list

    def add_condition(self, name: str, condition_list: List[str]):
        self._conditions[name] = condition_list

    def add_design(self, design: str):
        if len(design.strip()) == 0:
            self.log.error("Design can not be empty string")
        self.design = design

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

        # Create analysis folder if not exists
        make_path(self.config.deseq2_final_output)

    def _get_run_input(self, srr: str, meta: MetaParser):
        try:
            data = meta.data[srr]
            f = ""
            if self.method == "salmon":
                f = f"{data[meta.key_salmon][meta.key_output]}/quant.sf"
            elif self.method == "kallisto":
                f = f"{data[meta.key_kallisto][meta.key_output]}/abundance.tsv"
            elif self.method == "stringtie":
                f = f"{self.config.names.sra.folder(meta.sra)}" \
                    f"/{meta.key_string_tie}/t_data.ctab"
            elif self.method == "star":
                f = f"{data[meta.key_string_tie][meta.key_input]}"
            else:
                self.log.error(f"Method {self.method} is not supported for "
                               f"DESeq2 analysis yet.")

            if not exists_path(f):
                self.log.error(f"Input file not found : {f}")

            return f
        except KeyError:
            self.log.error(f"Input needed for DESeq2 analysis using "
                           f"{self.method} not found.")

    def _generate_config_file(self) -> str:
        files = {}
        all_runs = []
        run_types = []
        for sample in self.samples:
            files[sample] = {}
            meta = MetaParser(sample,
                              self.config.names.sra.folder(sample),
                              self.log)

            for rns in meta.get_runs():
                all_runs.append(rns)
                files[sample][rns] = self._get_run_input(rns, meta)
                run_types.append(int(meta.is_paired_end(rns)))

        mf = self._generate_meta_file(all_runs)

        # Following meta data is passed to the R scripts for the downstream
        # analysis

        data = {
            "input": files,
            "output": self.config.deseq2_final_output,
            "samples": self.samples,
            "paired": run_types,
            "conditions": self._conditions,
            "design": self.design,
            "method": self.method,
            "gtf": self.config.names.genome.gtf,
            "meta": mf,
            "contrast": self.contrast,
            "alpha": self.alpha,
            "threads": self.config.no_of_threads
        }
        jf = f"{self.config.deseq2_final_output}/{self.method}.config.json"
        with open(jf, "w", encoding='utf-8') as f:
            json.dump(data, f, indent=4)

        return jf

    def _generate_meta_file(self, runs) -> str:
        data = [runs]
        col_names = ["samples"]
        for condition in self._conditions:
            col_names.append(condition)
            data.append(self._conditions[condition])

        df_data = {}
        for c, d in zip(col_names, data):
            df_data[c] = d

        df = pd.DataFrame(df_data)
        meta_file = f"{self.config.deseq2_final_output}/{self.method}.meta.csv"
        df.to_csv(meta_file, index=False)
        self.log.info(f"Meta file saved in {meta_file}")
        return meta_file

    def _perform_analysis(self):
        f = self._generate_config_file()
        opts = [
            "Rscript",
            "rscripts/combined.R",
            "-c",
            f
        ]

        if subprocess.run(opts).returncode != 0:
            self.log.error("Something went wrong in running DESeq2 analysis")

    def run(self):
        self.config.log.info(f"DESeq2 analysis with {self.method} started")

        # Check setup
        self._check_setup()

        # analyse
        self._perform_analysis()

        self.config.log.info(f"DESeq2 analysis with {self.method} finished")
