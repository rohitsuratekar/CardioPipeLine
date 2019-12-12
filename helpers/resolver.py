#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# All path and filename resolver

import os
import shutil

from constants.genomics import *
from constants.system import *
from constants.tools import *


class BaseResolver:

    @staticmethod
    def exists(path: str) -> bool:
        """Checks if given path exists
        :param path: full path
        :return: True if path exists
        """
        return os.path.exists(path)

    def make(self, path: str):
        """ Makes dir structure of given path

        :param path: Path of the dictionary
        """
        if not self.exists(path):
            os.makedirs(path)

    def remove(self, path):
        """ Removes everything at this path (including sub-folders and files)
        :param path: Full path to delete
        """
        if not self.exists(path):
            return
        if os.path.isfile(path):
            os.remove(path)
        else:
            shutil.rmtree(path)

    @staticmethod
    def move(origin, destination):
        shutil.move(origin, destination)


class ToolResolver:
    def __init__(self):
        self.prefetch = f"{TOOL_SRA}/prefetch"
        self.fasterq_dump = f"{TOOL_SRA}/fasterq-dump"
        self.sortmerna = f"{TOOL_SORT_ME_RNA}/sortmerna"
        self.star = f"{TOOL_STAR}/STAR"
        self.salmon = f"{TOOL_SALMON}/salmon"
        self.salmon_decoy = TOOL_SALMON_DECOY
        self.mashmap = f"{TOOL_MASHMAP}/mashmap"
        self.bedtools = f"{TOOL_BEDTOOLS}/bedtools"
        self.kallisto = f"{TOOL_KALLISTO}/kallisto"
        self.unmerge = TOOL_UNMERGE


class FileNameResolver(BaseResolver):
    def __init__(self):
        self.sra_db = f"{DEX_DATA}/sra"
        self._ncbi_sra = f"{FOLDER_NCBI_PUBLIC}/sra"

    def star_bam_prefix(self, sra_id: str):
        return f"{self.sra_db}/{sra_id}/{sra_id}_"

    @staticmethod
    def sortmerna_align(index_path: str) -> str:
        return f"{index_path}/out/aligned.fastq"

    @staticmethod
    def sortmerna_out(index_path: str) -> str:
        return f"{index_path}/out/other.fastq"

    @staticmethod
    def sortmerna_log(index_path: str) -> str:
        return f"{index_path}/out/aligned.log"

    def filtered_read(self, no: int, sra_id: str, make=False) -> str:
        if make:
            self.make(f"{self.sra_db}/{sra_id}/filtered")
        return f"{self.sra_db}/{sra_id}/filtered/{sra_id}_{no}.filtered.fastq"

    def filtering_log(self, sra_id: str, make=False):
        if make:
            self.make(f"{self.sra_db}/{sra_id}/filtered")
        return f"{self.sra_db}/{sra_id}/filtered/aligned.log"

    def sra_meta(self, sra_id: str) -> str:
        """
        :param sra_id: SRA id
        :return: full path of metadata file of given SRA id
        """
        return f"{self.sra_db}/{sra_id}/{sra_id}_meta.json"

    def srr(self, sra_id: str, srr_id: str) -> str:
        """
        Generate name of the sra file of the given SRA id and its run id
        :param sra_id: SRA id (experiment id)
        :param srr_id: SRR id (run id)
        :return: full path of the run.sra file
        """
        return f"{self.sra_db}/{sra_id}/{srr_id}.sra"

    def sra_folder(self, sra_id):
        return f"{self.sra_db}/{sra_id}"

    def ncbi_srr(self, srr_id: str) -> str:
        """
        :param srr_id: SRR id
        :return: default path of the downloaded file in temporary NCBI
        public folder
        """
        return f"{self._ncbi_sra}/{srr_id}.sra"

    def salmon_output_folder(self, sra_id: str, make=True):
        if make:
            self.make(f"{self.sra_db}/{sra_id}/salmon")
        return f"{self.sra_db}/{sra_id}/salmon"

    @staticmethod
    def salmon_index_info(salmon_index_path: str):
        return f"{salmon_index_path}/info.json"

    @staticmethod
    def salmon_decoy_index(salmon_index_path: str) -> str:
        return f"{salmon_index_path}/gentrome.fa"

    @staticmethod
    def salmon_decoy_keys(salmon_index_path: str) -> str:
        return f"{salmon_index_path}/decoys.txt"

    @staticmethod
    def kallisto_index(index_path: str) -> str:
        return f"{index_path}/kallisto.idx"


class IndexResolver(BaseResolver):
    def __init__(self):
        self._index = f"{DEX_WORK}/Genome/index"

    @property
    def star(self):
        return f"{self._index}/star"

    @property
    def sortmerna(self):
        return f"{self._index}/sortmerna"

    @property
    def salmon(self):
        return f"{self._index}/salmon"

    @property
    def kallisto(self):
        return f"{self._index}/kallisto"


class PathUtil(BaseResolver):
    def __init__(self):
        self.data = DEX_DATA
        self.sra_accession = SRA_ACCESSION
        self.genome_fasta = ZEBRAFISH_GENOME_FASTA
        self.gtf_annotation = ZEBRAFISH_GTF_ANNOTATION
        self.transcript_fasta = ZEBRAFISG_TRANSCRIPT
        self.rrna_db = RRNA_DATABASE
        self._tools = None
        self._files = None
        self._index = None

    @property
    def tools(self) -> ToolResolver:
        if self._tools is None:
            self._tools = ToolResolver()
        return self._tools

    @property
    def files(self) -> FileNameResolver:
        if self._files is None:
            self._files = FileNameResolver()
        return self._files

    @property
    def index(self) -> IndexResolver:
        if self._index is None:
            self._index = IndexResolver()
        return self._index

    def sra_meta_exists(self, sra_id: str) -> bool:
        """
        Check if metadata file is present in the SRA folder

        This metadata file is generated by the CardioPipeLine package from
        the SRA_Annotation.tab file. It contains the relationship between
        given SRA id (experiment id) and the SRR id (run id). It also
        contains other metadata generated by this package as well as derived
        from the annotation file.
        :param sra_id: SRA id
        :return: returns True if metadata file exist for given SRA id
        """
        return self.exists(self.files.sra_meta(sra_id))

    def srr_exists(self, sra_id: str, srr_id: str) -> bool:
        """Checks if given sra file is available for given SRA id and SRR id
        :param sra_id: SRA id (experiment id)
        :param srr_id: SSR id (run id)
        :return:
        """
        return self.exists(self.files.srr(sra_id, srr_id))

    def fastq(self, sra_id) -> list:
        """
        :param sra_id: SRA id
        :return: List of all available fastq files in the folder of given
        SRA id
        """
        all_files = []
        for f in os.listdir(f"{self.files.sra_db}/{sra_id}"):
            if f.endswith(".fastq"):
                all_files.append(f"{self.files.sra_db}/{sra_id}/{f}")
        return all_files

    def fastq_dir(self, sra_id) -> str:
        """
        :param sra_id: SRA id
        :return: Path of fastq files of given SRA id
        """
        return f"{self.files.sra_db}/{sra_id}"

    def make_sra_path(self, sra_id: str):
        """Makes folder for given SRA id

        :param sra_id: SRA id
        """
        target = f"{self.files.sra_db}/{sra_id}"
        if not self.exists(target):
            os.makedirs(target)
