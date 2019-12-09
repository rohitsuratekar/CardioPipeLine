#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# All utils related to this project

from constants.system import *
from constants.tools import *
import os


class PathUtil:
    def __init__(self):
        self.data = DEX_DATA
        self.prefetch = f"{TOOL_SRA}/prefetch"
        self._ncbi_sra = f"{FOLDER_NCBI_PUBLIC}/sra"

    @property
    def sra(self):
        return f"{self.data}/sra"

    def sra_file(self, sra_id: str) -> str:
        return f"{self.sra}/{sra_id}/{sra_id}.sra"

    @staticmethod
    def exists(path: str) -> bool:
        return os.path.exists(path)

    def sra_exists(self, sra_id: str) -> bool:
        return self.exists(self.sra_file(sra_id))

    def ncbi_sra(self, sra_id: str):
        return f"{self._ncbi_sra}/{sra_id}.sra"
