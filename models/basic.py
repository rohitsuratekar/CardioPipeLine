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
        self.de_methods = []
        self._last_task = 16

    @property
    def data(self) -> dict:
        if self._data is None:
            with open(self.filename) as f:
                self._data = yaml.load(f, Loader=yaml.SafeLoader)
        return self._data

    @property
    def base(self) -> str:
        return self.data["base"]

    def _generate_tasks(self):
        self.tasks = []
        if isinstance(self.data['tasks'], int):
            if self.data['tasks'] == TASK_ALL:
                self.tasks = list(range(1, self._last_task + 1))
            else:
                self.tasks.append(self.data['tasks'])
        else:
            tsk = self.data['tasks'].split(',')
            correct_tasks = []
            for t in tsk:
                try:
                    correct_tasks.append(int(t))
                except ValueError as e:
                    r_range = str(t).strip().split("-")
                    if len(r_range) <= 1:
                        raise ValueError(f"'tasks' should be integer, list of "
                                         f"integers or range. e.g. tasks: 0 or "
                                         f"tasks: 1, 2, 3 ... or 2-5") from e
                    elif len(r_range) > 2:
                        raise ValueError(f"Range should be in "
                                         f"START_TASK-END_TASK format."
                                         f" e.g. 1-5, 11-15") from e
                    else:
                        try:
                            tmp = sorted([int(x) for x in r_range])
                            correct_tasks.extend(
                                list(range(tmp[0], tmp[1] + 1)))
                        except ValueError:
                            raise ValueError(
                                f"Range should have only integers, "
                                f"you have provided '{r_range}'")

            self.tasks = sorted(list(set(correct_tasks)))

    def _generate_methods(self):
        temp = []
        methods = ['star', 'salmon', 'stringtie', 'kallisto']
        for m in methods:
            if str(self.data['deseq2']['counts'][m]).lower() == "true":
                temp.append(m)
        if len(temp) == 0:
            raise ValueError("At lest one 'counts' method should be 'true' "
                             "for DESeq2 analysis. Please check your "
                             "'config.yaml' file.")
        self.de_methods = temp

    def is_paired(self, srr_id) -> bool:
        df = pd.read_csv(self.data['samples'])
        for _, row in df.iterrows():
            runs = row["run"]
            if runs.strip() == srr_id:
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

    def get_deseq2_inputs(self, method):
        condition_col = self.data['deseq2']['design_column']
        all_runs = []
        all_conditions = []
        df = pd.read_csv(self.data['samples'])
        for _, row in df.iterrows():
            all_runs.append(row["run"])
            all_conditions.append(row[condition_col])

        files = []
        for srr_id in all_runs:
            if method == "star":
                files.append(
                    f"deseq2/counts/{srr_id}/{srr_id}.star.counts")
            elif method == "stringite":
                files.append(
                    f"deseq2/counts/{srr_id}/{srr_id}.stringtie.counts")
            elif method == "salmon":
                files.append(
                    f"deseq2/counts/{srr_id}/{srr_id}.salmon.counts")
            elif method == "kallisto":
                files.append(
                    f"deseq2/counts/{srr_id}/{srr_id}.kallisto.counts")
            elif method == "stringtie":
                files.append(
                    f"deseq2/counts/{srr_id}/{srr_id}.stringtie.counts")
            else:
                raise ValueError(f"Method '{method}' is not supported for "
                                 f"DESeq2 analysis yet.")

        files = [f"{self.base}/{x}" for x in files]
        return files, all_runs, all_conditions

    def get_deseq2_outputs(self):
        self._generate_methods()
        files = []
        for m in self.de_methods:
            files.append(f"{self.base}/deseq2/analysis/{m}/analysis.log")
        return files

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
                    f"deseq2/counts/{srr_id}/{srr_id}.stringtie.counts")
                out.append(
                    f"deseq2/counts/{srr_id}/{srr_id}.stringtie_transcripts.counts")
            elif t == TASK_SALMON_INDEX:
                out.append("index/salmon/pos.bin")
            elif t == TASK_SALMON_MAPPING:
                out.append(f"mappings/salmon/{srr_id}/quant.sf")
            elif t == TASK_KALLISTO_INDEX:
                out.append("index/kallisto/kallisto.idx")
            elif t == TASK_KALLISTO_MAPPING:
                out.append(f"mappings/kallisto/{srr_id}/abundance.tsv")
            elif t == TASK_COUNT_STAR:
                out.append(f"deseq2/counts/{srr_id}/{srr_id}.star.counts")
            elif t == TASK_COUNT_SALMON:
                out.append(f"deseq2/counts/{srr_id}/{srr_id}.salmon.counts")
            elif t == TASK_COUNT_KALLISTO:
                out.append(f"deseq2/counts/{srr_id}/{srr_id}.kallisto.counts")
            elif t < self._last_task + 1:
                pass
            else:
                raise KeyError(f"Invalid task ID: '{t}'. Check README.md "
                               f"file to know what are the task IDs available")

        return [f"{self.base}/{x}" for x in out]


def run(self):
    d = self.fastq("SRR6039671")
    print(d)
