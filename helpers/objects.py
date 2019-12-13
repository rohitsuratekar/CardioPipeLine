#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# All objects will go here

import os
import subprocess
import sys

from helpers.logging import Log
from helpers.utils import make_path


class Tool:
    def __init__(self, path: str, binary: str, log: Log):
        self._path = path
        self._binary = binary
        self._output_folder = None
        self._log = log
        self._index = None

        self._log.debug(f"{self.log_header} Tool is initialized")

    @property
    def path(self) -> str:
        return f"{self._path}/{self._binary}"

    @property
    def log_header(self) -> str:
        return f"[ {self._binary} ]"

    @property
    def output_folder(self):
        if self._output_folder is None:
            self._output_folder = self._default_path("Output Folder")
        return self._output_folder

    @property
    def index(self):
        if self._index is None:
            self._index = self._default_path("Index")
        return self._index

    @index.setter
    def index(self, path: str):
        self._index = path

    @output_folder.setter
    def output_folder(self, path: str):
        self._output_folder = path

    def _default_path(self, use: str):
        # If folder is not defined, use script folder
        path = os.path.dirname(os.path.realpath(sys.argv[0]))
        path = f"{path}/{self._binary.strip().lower()}"
        self._log.warn(f"{self.log_header} Default path for {use} is not "
                       f"defined and hence using base script folder "
                       f"{path}")
        return path

    def _validate(self, options: list):
        for opt in options:
            if type(opt) != str:
                self._log.error(f"All the options provided should be String. "
                                f"Please wrap option '{opt}' in string.")

    def run(self, options: list, success: str = None, error: str = None,
            exception=None):
        self._validate(options)
        self._log.debug(
            f"{self.log_header} {self._binary} {' '.join(options)}")
        opt = [self.path]  # Add binary
        opt.extend(options)  # All all options
        if subprocess.run(opt).returncode == 0:
            if success is not None:
                self._log.info(success)
            else:
                self._log.info(f"{self.log_header} Success")
        else:
            if error is not None:
                self._log.error(error, exception=exception)
            else:
                self._log.error(f"{self.log_header} Something went wrong")


class Genome:
    def __init__(self, data: dict):
        self.fasta = data["fasta"]
        self.gtf = data["gtf_annotation"]
        self.transcript = data["transcript"]
        self.rrna_db = data["rrna_database"]


class SRA:
    def __init__(self, data: dict, log: Log):
        self.sra_db = data["path"]
        self.accession = data["accession_db"]
        self.ncbi_public = data["ncbi_public"]
        self._log = log

    def folder(self, sra_id, make=False) -> str:
        path = f"{self.sra_db}/{sra_id}"
        if make:
            make_path(path)
        return path

    def meta(self, sra_id: str) -> str:
        return f"{self.folder(sra_id)}/{sra_id}_meta.json"


class NameResolver:
    def __init__(self, data: dict, log: Log):
        self.data = data
        self._log = log
        self._sra = None
        self._genome = None

    @property
    def genome(self) -> Genome:
        if self._genome is None:
            self._genome = Genome(self.data["genome"])
        return self._genome

    @property
    def sra(self) -> SRA:
        if self._sra is None:
            self._sra = SRA(self.data["sra"], self._log)
        return self._sra
