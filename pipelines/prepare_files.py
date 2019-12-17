#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# Pipeline to prepare files for specific analysis

import os

from helpers import ConfigParser, move_path, exists_path, MetaParser
from pipelines._base_pipe import PipeLine


class PrepareRNAseq(PipeLine):
    def __init__(self, config: ConfigParser, meta: MetaParser):
        super().__init__()
        self.sra = meta.sra
        self.config = config
        self.meta = meta

        self.log.info("PrepareRNAseq pipeline initialized")

    def download_run(self, srr: str):
        opts = [
            "--verbose",  # Verbose to see logs
            srr.strip()  # Run ID for downloading
        ]
        self.log.info("Downloading started")
        self.config.tools.prefetch.run(opts,
                                       success=f"{srr} run is downloaded "
                                               f"in {self.config.names.sra.ncbi_public}")

        # Move to the SRA folder for further analysis
        origin = f"{self.config.names.sra.ncbi_public}/sra/{srr}.sra"
        destination = self.config.names.sra.raw_data(self.sra, srr)
        move_path(origin, destination)
        self.log.info(f"SRA file for {srr} downloaded and moved to default "
                      f"SRA folder {destination}")

    def check_and_download_raw_data(self) -> bool:
        runs = self.meta.get_runs()
        if len(runs) == 0:
            self.log.info(f"No runs found for {self.sra}. Trying to extract "
                          f"them from the SRA accession database")
            self.meta.extract_runs(self.config.names.sra.accession)
            runs = self.meta.get_runs()

        for r in runs:
            if not exists_path(self.config.names.sra.raw_data(self.sra, r)):
                self.log.info(f"Run {r} runs not found in the SRA "
                              f"folder."
                              f" Downloading them with prefetch")
                self.download_run(r)

        self.log.info(f"All the runs in {self.sra} extracted from "
                      f"{self.config.names.sra.folder(self.sra)}")
        return True

    def convert_to_fastq(self, ssr: str):
        opts = [
            self.config.names.sra.raw_data(self.sra, ssr),  # Path to
            # input SRA file
            "--split-files",  # Splits files according to Reads
            "--progress",  # Show progress
            "--skip-technical",  # Skip technical reads
            "--threads",  # Number of threads
            str(self.config.no_of_threads),
            # If you are unsure, use 6 (default)
            "--outdir",  # Output folder
            self.config.names.sra.folder(self.sra)
        ]

        self.config.tools.fasterq_dump.run(opts,
                                           success=f"{ssr} is converted "
                                                   f"into the fastq files")

        # Add converted fastq files to the metadata
        files = []
        folder = self.config.names.sra.folder(self.sra)
        for file in os.listdir(folder):
            if os.path.isfile(f"{folder}/{file}"):
                if str(file).endswith(".fastq") and str(
                        file).startswith(ssr):
                    files.append(f"{folder}/{file}")
        data = self.meta.data
        data[ssr][self.meta.key_fastq] = sorted(files)
        data[ssr][self.meta.key_pair_end] = len(files) == 2
        self.meta.append(data)

    def check_if_fastq_exists(self, ssr) -> bool:
        fastq = self.meta.get_fastq(ssr)
        if len(fastq) == 0:
            self.log.info(f"Fastq files not found for {self.sra}")
            return False

        self.log.info(f"Fastq files collected from "
                      f"{self.config.names.sra.folder(self.sra)}")
        return True

    def check_metadata(self):
        if not exists_path(self.meta.path):
            self.log.info("There is no metadata file. Generating one...")
            self.meta.extract_runs(self.config.names.sra.accession)
            self.log.info(f"Metadata file for {self.sra} generated")

    def run(self):
        # Check if meta data exist and if not, create one
        self.check_metadata()

        # Check if filtered files exits and then return

        for srr in self.meta.get_runs():

            filtered = self.meta.get_filtered(srr)
            filtered = [exists_path(x) for x in filtered]
            if all(filtered) and len(filtered) != 0:
                self.log.info(f"Filtered files exists for {srr}")
                continue

            # Check if fastq file exists
            if not self.check_if_fastq_exists(srr):
                # IF no fastq files are there, download them
                self.check_and_download_raw_data()
                # and convert them again
                self.convert_to_fastq(srr)

        self.log.info(f"All fastq files for {self.sra} are prepared for "
                      f"downstream analysis")
        self.log.info("Exiting PrepareRNAseq pipeline")
