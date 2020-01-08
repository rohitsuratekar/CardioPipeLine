#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# DESeq2 pipeline


import subprocess
from typing import Union

from helpers import ConfigParser, MetaParser, make_path, exists_path
from pipelines._base_pipe import PipeLine


class DESeq2(PipeLine):
    def __init__(self, config: ConfigParser, meta: MetaParser, method: str):
        super().__init__()
        self.config = config
        self.meta = meta
        self.sra = meta.sra
        self.method = method.strip().lower()
        self.base_path = f"{self.config.names.sra.folder(self.sra)}/deseq2"
        self.log.info(f"DESeq2 pipeline started for method : {method}")

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

        new_data = {
            "counts": {
                self.method: self.counts_file(srr)
            }
        }

        data["counts"] = {**data["counts"], **new_data["counts"]}
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

    def _prepare(self):
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
        self._prepare()
