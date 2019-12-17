#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# rRNA filtering pipeline

import os
import subprocess
import sys

from helpers import (ConfigParser, MetaParser,
                     delete_path, exists_path, move_path, make_path)
from pipelines._base_pipe import PipeLine


class RRNAFiltering(PipeLine):
    def __init__(self, config: ConfigParser, meta: MetaParser):
        super().__init__()
        self.sra_id = meta.sra
        self.config = config
        self.meta = meta

        self.log.info("rRNAFiltering pipeline initialized")

    def _check_if_exists(self, srr):
        filtered = self.meta.get_filtered(srr)
        for f in filtered:
            if not exists_path(f):
                return False
        self.log.info(f"Filtered reads for {srr} already exists")
        return True

    def perform_filtering(self):
        self.log.info("Attempting to generate index for SortMeRNA")
        for srr in self.meta.get_runs():

            if self._check_if_exists(srr):
                # If filtered files already exists, skip the loop
                continue

            opts = []
            db_path = self.config.names.genome.rrna_db
            for file in os.listdir(db_path):
                # TODO: Proper handling from the config file
                if str(file).endswith(".fasta"):
                    if str(file).startswith("silva-euk"):
                        opts.append("-ref")
                        opts.append(f"{db_path}/{file}")

            for file in self.meta.get_fastq(srr):
                opts.append("-reads")  # For input reads
                opts.append(file)

            # Add other options if "-paired_in" is used, both reads will be
            # merged, hence don't use it
            opts.extend([
                "-workdir",  # Working folder for storing index and output
                self.config.tools.sortmerna.index,
                "-fastx",  # Output reads into fastq format
                "-other",
                # Output the Other.fastq file (which has NON-rRNA reads)
                "-a",  # Number of threads
                str(self.config.no_of_threads),  # If not available, use 1
                "-v",  # Produce verbose output
                "-num_alignments",  # Only first alignment passing of E value
                str(1),
                # This will be fastest choice when only filtering is needed,
                "-paired_in"
                # if either of read match, put both reads into the
                # 'aligned.fasta' file.
            ])

            # Delete the kvdb folder (needed for SortMeRNA)
            delete_path(f"{self.config.tools.sortmerna.index}/kvdb")

            self.config.tools.sortmerna.run(opts,
                                            success="rRNA filtering with "
                                                    "SortMeRNA is successful")

            self.post_processing(srr)

    def post_processing(self, srr: str):
        out_file = self.config.tools.sortmerna.get_extra("out")
        rna_file = self.config.tools.sortmerna.get_extra("align")
        log_file = self.config.tools.sortmerna.get_extra("log")

        make_path(f"{self.config.names.sra.folder(self.sra_id)}/"
                  f"{self.config.filtered_id}")

        # Sanity check
        if not exists_path(out_file):
            self.log.error(
                f"Output file for SortMeRNA are not found at {out_file}")

        base_folder = os.path.dirname(sys.argv[0])
        unmerg_script = f"{base_folder}/external/unmerge-paired-reads.sh"

        # Check if binary has proper access to execute the command
        if not os.access(unmerg_script, os.X_OK):
            self._log.error(f"Script {unmerg_script} do not "
                            f"have permission to execute.",
                            exception=PermissionError)

        opts = [
            unmerg_script,
            out_file,
            self.config.names.sra.filtered_srr(self.sra_id, srr, 1),
            self.config.names.sra.filtered_srr(self.sra_id, srr, 2)
        ]

        if subprocess.run(opts).returncode == 0:
            self.log.info("Filtered files successfully un-merged")
        else:
            self.log.error("Something went wrong in un-merging filtered files")

        # Move log file to filtered folder
        move_path(log_file, self.config.names.sra.filtered_log(self.sra_id,
                                                               srr))

        # Delete filtered files
        delete_path(rna_file)
        delete_path(out_file)

        self.log.info(f"rRNA filtering post processing done for {srr}")

        # Add data to the metadata file

        data = self.meta.data
        data[srr][self.meta.key_filtered] = [
            self.config.names.sra.filtered_srr(self.sra_id, srr, 1),
            self.config.names.sra.filtered_srr(self.sra_id, srr, 2)
        ]
        self.meta.append(data)

    def run(self):
        if not self.config.rrna_filtering:
            self.log.info("rRNA filtering is switched OFF in the config.json "
                          "file. Please use 'filter_rrna=1' and start again.")
            self.log.info("Exiting rRNAFiltering pipeline")
            return

        self.perform_filtering()
        self.log.info("Exiting rRNAFiltering pipeline")
