#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# All parsers will go here

import json
import logging
import os

from helpers.logging import Log
from helpers.objects import Tool, NameResolver
from helpers.utils import make_path, exists_path


class ToolParser:
    def __init__(self, data: dict, log: Log, index_path: str):
        self.data = data
        self.index_path = index_path
        self._sra_tool = None
        self._fasterq_dump = None
        self._star = None
        self._salmon = None
        self._kallisto = None
        self._mashmap = None
        self._bedtools = None
        self._sortmerna = None
        self._log = log

    def _make_tool(self, key: str, binary: str) -> Tool:
        path = self.data[key]
        # Check if binary exists in given path
        if not os.path.exists(f"{path}/{binary}"):
            self._log.error(
                f"'{binary}' does not exist in path {path}."
                f" Please check path of your tool in 'config.json' file",
                exception=FileNotFoundError)

        # Check if binary has proper access to execute the command
        if not os.access(f"{path}/{binary}", os.X_OK):
            self._log.error(f"'{binary}' at {path} do not "
                            f"have permission to execute.",
                            exception=PermissionError)

        return Tool(path, binary, self._log)

    @property
    def prefetch(self) -> Tool:
        if self._sra_tool is None:
            self._sra_tool = self._make_tool("sratools", "prefetch")
        return self._sra_tool

    @property
    def fasterq_dump(self) -> Tool:
        if self._fasterq_dump is None:
            self._fasterq_dump = self._make_tool("sratools", "fasterq-dump")
        return self._fasterq_dump

    @property
    def star(self) -> Tool:
        if self._star is None:
            self._star = self._make_tool("star", "STAR")
            self._star.index = f"{self.index_path}/star"
        return self._star

    @property
    def salmon(self) -> Tool:
        if self._salmon is None:
            self._salmon = self._make_tool("salmon", "salmon")
            self._salmon.index = f"{self.index_path}/salmon"
        return self._salmon

    @property
    def kallisto(self) -> Tool:
        if self._kallisto is None:
            self._kallisto = self._make_tool("kallisto", "kallisto")
            self._kallisto.index = f"{self.index_path}/kallisto"
        return self._kallisto

    @property
    def mashmap(self) -> Tool:
        if self._mashmap is None:
            self._mashmap = self._make_tool("mashmap", "mashmap")
        return self._mashmap

    @property
    def bedtools(self) -> Tool:
        if self._bedtools is None:
            self._bedtools = self._make_tool("bedtools", "bedtools")
        return self._bedtools

    @property
    def sortmerna(self) -> Tool:
        if self._sortmerna is None:
            self._sortmerna = self._make_tool("sortmerna", "sortmerna")
            self._sortmerna.index = f"{self.index_path}/sortmerna"
            self._sortmerna.output_folder = f"{self._sortmerna.index}/out"
        return self._sortmerna


class MetaParser:
    def __init__(self, sra: str, folder: str, log: Log):
        self.sra = sra
        self.folder = folder
        self.name = f"{sra}_meta.json"
        self._data = {}  # type: dict
        self._log = log

        self._key_runs = "Runs"
        self._key_fastq = "fastq"
        make_path(folder)
        # Check if file exists
        if exists_path(f"{folder}/{self.name}"):
            log.info(f"{sra} metadata exists, loading data.")
            with open(f"{folder}/{self.name}") as f:
                self._data = json.load(f)

    @property
    def data(self) -> dict:
        return self._data

    @property
    def path(self):
        return f"{self.folder}/{self.name}"

    def save(self):
        with open(self.path, "w", encoding='utf-8') as f:
            json.dump(self.data, f, indent=4)
        self._log.info(f"{self.sra} Metadata updated and saved")

    def append(self, data: dict):
        self._data = {**self.data, **data}
        self.save()

    def get_runs(self) -> list:
        if self._key_runs in self.data.keys():
            return self.data[self._key_runs]
        else:
            return []

    def get_raw(self, srr_id: str) -> list:
        if srr_id not in self.data.keys():
            self._log.error(f"Run '{srr_id}' not found in {self.sra}",
                            exception=FileNotFoundError)
        return self.data[srr_id][self._key_fastq]


class ConfigParser:
    def __init__(self, config: dict):
        self.data = config
        self._tools = None
        self._names = None
        self._log = None

    @property
    def log(self) -> Log:
        debug = -1
        log_format = self.data["setup"]["log_format"]
        if not self.data["setup"]["show_debug_log"]:
            debug = logging.DEBUG

        if len(str(log_format).strip()) == 0:
            log_format = "%(message)s"
        if self._log is None:
            log = Log(show_log=self.data["setup"]["show_log"],
                      add_to_file=self.data["setup"]["save_log"],
                      filename=self.data["setup"]["log_file"],
                      logging_format=log_format,
                      min_log_level=debug)
            self._log = log
            self._log.info("Logging is initialized")
        return self._log

    @property
    def tools(self) -> ToolParser:
        if self._tools is None:
            self._tools = ToolParser(self.data["tools"], self.log,
                                     self.data["index_path"])
        return self._tools

    @property
    def names(self) -> NameResolver:
        if self._names is None:
            self._names = NameResolver(self.data, self.log)
        return self._names
