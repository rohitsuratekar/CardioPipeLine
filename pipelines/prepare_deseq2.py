#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# DESeq2 pipeline


import subprocess
from typing import Union

import pandas as pd

from helpers import ConfigParser, MetaParser, make_path, exists_path
from helpers.constants import CM_GENE_ID
from pipelines._base_pipe import PipeLine


class CountMatrix(PipeLine):
    def __init__(self, config: ConfigParser, meta: MetaParser, method: str):
        super().__init__()
        self.config = config
        self.meta = meta
        self.sra = meta.sra
        self.method = method.strip().lower()
        self.base_path = f"{self.config.names.sra.folder(self.sra)}/deseq2"
        self.log.info(f"CountMatrix pipeline started for method : {method}")

    def _get_value(self, *args) -> Union[str, None]:
        # Checks if given path is present in the nested json of metadata,
        # if present, returns it else it will return None
        try:
            data = self.meta.data
            for k in args:
                data = data[k]
            return data
        except KeyError:
            return None

    def counts_file(self, srr: str) -> str:
        return f"{self.base_path}/{srr}_countmatrix.{self.method}.csv"

    def _run_script(self, name: str, arguments: list = None):
        opt = ["Rscript", f"rscripts/{name}"]
        if arguments is not None:
            opt.extend(arguments)
        if subprocess.run(opt).returncode == 0:
            self.log.info(f"R Script {name} finished")
        else:
            self.log.error(f"Problem occurred while running script {name}")

    def _save_meta(self, srr: str):
        data = {"counts": {}}
        if self.meta.key_deseq2 in self.meta.data[srr].keys():
            data = self.meta.data[srr][self.meta.key_deseq2]

        c = self.meta.key_counts
        new_data = {
            c: {
                self.method: self.counts_file(srr)
            }
        }
        data[c] = {**data[c], **new_data[c]}
        nd = self.meta.data
        nd[srr][self.meta.key_deseq2] = data
        self.meta.append(nd)

    def _prepare_star(self):
        # Prepare count-matrix
        for srr in self.meta.get_runs():

            if exists_path(self.counts_file(srr)):
                self.log.info(f"Countmatrix for {srr} is already generated "
                              f"at {self.counts_file(srr)}")
                continue

            bam_file = self._get_value(srr, self.meta.key_string_tie,
                                       self.meta.key_input)
            if bam_file is None:
                self.log.error(f"BAM file not found in the {srr} metadata "
                               f"file. Please run StringTie/STAR pipeline "
                               f"before running DESeq2 with 'star'")
            opts = [
                "-f",
                bam_file,
                "-g",
                self.config.names.genome.gtf,
                "-s",
                srr,
                "-o",
                self.counts_file(srr),
                "-p",
                str(int(self.meta.is_paired_end(srr))),
                "-t",
                str(self.config.no_of_threads)
            ]
            self._run_script("countmatrix.R", opts)

            self.log.info(f"Countmatrix generated at {self.counts_file(srr)}")
            self._save_meta(srr)

    def _tximport_input(self, srr: str):
        if self.method == "salmon":
            return f"{self._get_value(srr, self.meta.key_salmon, self.meta.key_output)}/quant.sf"
        elif self.method == "kallisto":
            return f"{self._get_value(srr, self.meta.key_kallisto, self.meta.key_output)}/abundance.tsv"
        elif self.method == "stringtie":
            return f"{self.config.names.sra.folder(self.sra)}/" \
                   f"{self.meta.key_string_tie}/t_data.ctab"
        else:
            self.log.error(f"Input method {self.method} is not supported yet")

    def _prepare_tximport(self):

        for srr in self.meta.get_runs():
            if exists_path(self.counts_file(srr)):
                self.log.info(f"Countmatrix for {srr} is already generated "
                              f"at {self.counts_file(srr)}")
                continue

            in_file = self._tximport_input(srr)

            if in_file is None:
                self.log.error(f"Input file from {self.method} analysis not "
                               f"found for {srr}. Please run {self.method} "
                               f"Pipeline before DESeq2")
            opts = [
                "-m",
                self.method,
                "-f",
                in_file,
                "-s",
                srr,
                "-g",
                self.config.names.genome.gtf,
                "-o",
                self.counts_file(srr)
            ]
            self._run_script("tximport.R", opts)
            self.log.info(f"Countmatrix generated at {self.counts_file(srr)}")
            self._save_meta(srr)

    def _prepare_cm(self):
        if self.method == "star":
            self._prepare_star()
        elif self.method in ["salmon", "kallisto", "stringtie"]:
            self._prepare_tximport()
        else:
            self.log.error(f"Method {self.method} is not supported by DESeq2 "
                           f"pipeline yet.")

    def run(self):
        # Create DESeq2 folder for storing the output files
        make_path(self.base_path)
        # Prepare the input data for the analysis
        self._prepare_cm()

        self.log.info(f"CountMatrix pipeline finished for {self.sra}")


class DESeqMeta(PipeLine):

    def __init__(self, config: ConfigParser, method: str):
        super().__init__()
        self.config = config
        self.method = method.strip().lower()
        self._samples = []
        self._conditions = {}
        self.log.info(f"GenerateMeta pipeline started for method : {method}")

    def samples(self, samples: list):
        # Sanity check to avoid spaces and empty ids
        self._samples = [x.strip() for x in samples if len(x.strip()) > 0]
        self.log.info(f"Samples added to the analysis: {self._samples}")
        return self

    def add_conditions(self, name: str, details: list):
        self._conditions[name.strip()] = details
        self.log.info(f"Condition {name} added to the analysis")
        return self

    def get_meta(self, sra: str) -> MetaParser:
        return MetaParser(sra, self.config.names.sra.folder(sra),
                          self.config.log)

    def _check_matrix(self) -> dict:
        files = {}
        for sample in self._samples:
            meta = self.get_meta(sample)
            for srr in meta.get_runs():
                try:
                    f = meta.data[srr][meta.key_deseq2][meta.key_counts][
                        self.method]
                    if not exists_path(f):
                        self.log.error(
                            f"Count matrix generated by {self.method}"
                            f" not found for {srr} of {sample}. Please "
                            f"run CountMatrix pipeline before "
                            f"generating meta data.")
                    files[sample] = {srr: f}
                except KeyError:
                    self.log.error(f"Count matrix details are not present in "
                                   f"the metadata file. Please run "
                                   f"CountMatrix pipeline before generating "
                                   f"meta data for DESeq2.")
        return files

    def _combine_files(self, file_list: list) -> str:
        dfs = []
        for f in file_list:
            dfs.append(pd.read_csv(f))

        dfs = [x.set_index(CM_GENE_ID) for x in dfs]
        dfs = pd.concat(dfs, axis=1, join="outer", sort=False)
        dfs = dfs.sort_index().reset_index()
        count_file = f"{self.config.deseq2_final_output}/counts.{self.method}.csv"
        dfs.to_csv(count_file, index=False)
        self.log.info(f"Final count file saved in {count_file}")
        return count_file

    def _check_setup(self):
        samples = len(self._samples)
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

    def _create_meta_file(self) -> str:
        data = [self._samples]
        col_names = ["samples"]
        for condition in self._conditions:
            col_names.append(condition)
            data.append(self._conditions[condition])

        df_data = {}
        for c, d in zip(col_names, data):
            df_data[c] = d

        df = pd.DataFrame(df_data)
        meta_file = f"{self.config.deseq2_final_output}/meta.{self.method}.csv"
        df.to_csv(meta_file, index=False)
        self.log.info(f"Meta file saved in {meta_file}")
        return meta_file

    def run(self):
        # Check analysis setup
        self._check_setup()

        # Create folder if does not exists
        make_path(self.config.deseq2_final_output)

        # Check if all the samples have count matrix and get their details
        samples = self._check_matrix()
        file_list = []
        rns = []
        for s in samples.values():
            for f in s:
                rns.append(f)
                file_list.append(s[f])

        # Merge the counts
        self.log.info(f"Merging runs from : {rns}")
        self._combine_files(file_list)

        # Create meta file
        self._create_meta_file()

        self.log.info("Exiting DESeqMeta pipeline.")
