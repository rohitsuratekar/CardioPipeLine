"""
CardioPipeLine (c) 2020

Author: Rohit Suratekar, ZDG Lab, IIMCB

This deals with all rules related STAR analysis
"""

import yaml


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
        if isinstance(self.data['tasks'], int):
            if self.data['tasks'] == 0:
                self.tasks = list(range(1, 15))
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

    def output_files(self, srr_id: str, is_paired: bool):
        self._generate_tasks()
        out = []
        for t in self.tasks:
            if t == 1:
                # Download SRA files
                out.append(f"sra/{srr_id}.sra")
            elif t == 2:
                # Convert to fastq
                if is_paired:
                    out.append(f"fastq/{srr_id}.sra_1.fastq")
                    out.append(f"fastq/{srr_id}.sra_2.fastq")
                else:
                    out.append(f"fastq/{srr_id}.sra.fastq")
            elif t == 3:
                if is_paired:
                    out.append(f"filtered/{srr_id}.sra.filtered_1.fastq")
                    out.append(f"filtered/{srr_id}.sra.filtered_2.fastq")
                else:
                    out.append(f"filtered/{srr_id}.sra.filtered.fastq")
            elif t == 4:
                out.append("index/star/SAindex")
            # elif t == 5:
            #     out.append(f"bams/{srr_id}.sra.Aligned.sortedByCoord.out.bam")
            elif t < 15:
                pass
            else:
                raise KeyError(f"Invalid task ID: '{t}'. Check README.md "
                               f"file to know what are the task IDs available")

        return [f"{self.base}/{x}" for x in out]

    def run(self):
        d = self.output_files("SRR4949", True)
        print(d)
