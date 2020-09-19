"""
CardioPipeLine (c) 2020

Author: Rohit Suratekar, ZDG Lab, IIMCB

This deals with PipeLine class which handles all the output file name
generations.
"""

import yaml
import pandas as pd
from models.constants import *


class PipeLine:
    def __init__(self, config_file: str):
        self.filename = config_file
        self._data = None
        self.tasks = []

    @property
    def data(self) -> dict:
        if self._data is None:
            with open(self.filename) as f:
                self._data = yaml.load(f, Loader=yaml.SafeLoader)
        return self._data

    @property
    def base(self) -> str:
        return self.data["base"]

    def _add_mandatory(self):
        # STAR mapping and STAR alignment
        if 5 in self.tasks and 4 not in self.tasks:
            self.tasks.append(4)

    def _generate_tasks(self):
        self.tasks = []
        if isinstance(self.data['tasks'], int):
            if self.data['tasks'] == TASK_ALL:
                self.tasks = list(range(1, 16))
            else:
                self.tasks.append(self.data['tasks'])
        else:
            tsk = self.data['tasks'].split(',')
            try:
                tsk = [int(x) for x in tsk]
            except ValueError as e:
                raise ValueError(f"'tasks' should be integer or list of "
                                 f"integers. e.g. tasks: 0 or tasks: 1, 2, "
                                 f"3 ...") from e
            self.tasks = tsk

        self._add_mandatory()

    def is_paired(self, srr_id) -> bool:
        df = pd.read_csv(self.data['samples'])
        for _, row in df.iterrows():
            run = row["run"]
            if run.strip() == srr_id:
                return row["is_paired"]
        raise KeyError(f"{srr_id} not found in 'samples.csv'")

    def fastq(self, srr_id):
        self._generate_tasks()
        if TASK_FILTER_RRNA in self.tasks:
            if self.is_paired(srr_id):
                return [
                    f"{self.base}/filtered/{srr_id}.sra.filtered_1.fastq",
                    f"{self.base}/filtered/{srr_id}.sra.filtered_2.fastq",
                ]
            else:
                return f"{self.base}/filtered/{srr_id}.sra.filtered.fastq"
        else:
            if self.is_paired(srr_id):
                return [
                    f"{self.base}/fastq/{srr_id}.sra_1.fastq",
                    f"{self.base}/fastq/{srr_id}.sra_2.fastq"
                ]
            else:
                return f"{self.base}/fastq/{srr_id}.sra.fastq"

    def output_files(self, srr_id: str, is_paired: bool):
        self._generate_tasks()
        out = []
        for t in self.tasks:
            if t == TASK_DOWNLOAD:
                # Download SRA files
                out.append(f"sra/{srr_id}.sra")
            elif t == TASK_CONVERT_FASTQ:
                # Convert to fastq
                if is_paired:
                    out.append(f"fastq/{srr_id}.sra_1.fastq")
                    out.append(f"fastq/{srr_id}.sra_2.fastq")
                else:
                    out.append(f"fastq/{srr_id}.sra.fastq")
            elif t == TASK_FILTER_RRNA:
                if is_paired:
                    out.append(f"filtered/{srr_id}.sra.filtered_1.fastq")
                    out.append(f"filtered/{srr_id}.sra.filtered_2.fastq")
                else:
                    out.append(f"filtered/{srr_id}.sra.filtered.fastq")
            elif t == TASK_STAR_INDEX:
                out.append("index/star/SAindex")
            elif t == TASK_STAR_MAPPING:
                out.append(f"bams/{srr_id}.sra.Aligned.sortedByCoord.out.bam")
            elif t == TASK_STRINGTIE:
                out.append(
                    f"mappings/stringtie/{srr_id}/{srr_id}_gene_expression.tsv")
                out.append(
                    f"mappings/stringtie/{srr_id}/{srr_id}_assembled.gtf")
            elif t == TASK_COUNT_STRINGTIE:
                out.append(
                    f"mappings/stringtie/{srr_id}/{srr_id}.gene_count_matrix.csv")
                out.append(
                    f"mappings/stringtie/{srr_id}/{srr_id}.transcript_count_matrix.csv")
            elif t == TASK_SALMON_INDEX:
                out.append("index/salmon/pos.bin")
            elif t == TASK_SALMON_MAPPING:
                out.append(f"mappings/salmon/{srr_id}/quant.sf")
            elif t == TASK_KALLISTO_INDEX:
                out.append("index/kallisto/kallisto.idx")
            elif t == TASK_KALLISTO_MAPPING:
                out.append(f"mappings/kallisto/{srr_id}/abundance.tsv")
            elif t == TASK_DESEQ:
                # TODO
                out.append("test.txt")
            elif t < 16:
                pass
            else:
                raise KeyError(f"Invalid task ID: '{t}'. Check README.md "
                               f"file to know what are the task IDs available")

        return [f"{self.base}/{x}" for x in out]

    def run(self):
        d = self.fastq("SRR6039671")
        print(d)
