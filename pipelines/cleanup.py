#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# Cleaning up files


from helpers import ConfigParser, MetaParser, delete_path
from pipelines._base_pipe import PipeLine


class CleanUp(PipeLine):
    def __init__(self, config: ConfigParser, meta: MetaParser):
        super().__init__()
        self.config = config
        self.meta = meta
        self.sra = meta.sra

        self.log.info("CleanUp pipeline initialized")

    def clean_up(self, srr: str):
        # Check if fastaq or filtered files are available
        fastq = self.meta.get_fastq(srr)
        filtered = self.meta.get_filtered(srr)

        # If filtered files are present, remove respective fastq files
        if len(filtered) > 0:
            for f in fastq:
                delete_path(f)
                self.log.info(f"Deleted fastq file : {f}")

            data = self.meta.data
            if self.meta.key_fastq in data.keys():
                del data[srr][self.meta.key_fastq]
                self.meta.append(data)

    def run(self):
        for srr in self.meta.get_runs():
            self.clean_up(srr)

        self.log.info("Exiting CleanUp pipeline")
