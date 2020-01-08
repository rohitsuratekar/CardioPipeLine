#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# STAR pipeline

import os
import sys

from helpers import ConfigParser, MetaParser, make_path, delete_path
from pipelines._base_pipe import PipeLine


class Star(PipeLine):
    def __init__(self, config: ConfigParser, meta: MetaParser,
                 force_index: bool = False):
        super().__init__()
        self.config = config
        self.meta = meta
        self.sra = meta.sra
        self.force_index = force_index

        output_path = f"{self.config.names.sra.folder(self.sra)}/star"
        make_path(output_path)
        self.config.tools.star.output_folder = output_path
        self.log.info(f"Star output folder set to {output_path}")

        self.log.info("STAR pipeline initialized")

    def index_star(self):
        # Check if index is already present

        index_path = self.config.tools.star.index
        make_path(index_path)

        if len(os.listdir(index_path)) > 0 and not self.force_index:
            self.log.info("STAR index folder is non-empty. Assuming index "
                          "already exists.")
            self.log.info("If you want to re-create index, please pass "
                          "'force_index=True' argument")
            return

        self.log.info("STAR indexing started")

        opts = [
            "--runMode",  # Define run mode
            "genomeGenerate",  # Mode to generate index
            "--runThreadN",  # Number of threads for running this program
            str(self.config.no_of_threads),  # If not known keep it as 1
            "--genomeDir",  # Folder where index should be made
            index_path,
            "--genomeFastaFiles",  # Fasta file of full genome
            self.config.names.genome.fasta,  # e.g.
            # Danio_rerio.GRCz11.dna.primary_assembly.fa
            "--sjdbGTFfile",  # GTF Annotation file
            self.config.names.genome.gtf  # e.g. Danio_rerio.GRCz11.98.chr.gtf
        ]

        self.config.tools.star.run(opts, success=f"STAR index generated "
                                                 f"successfully in {index_path}")

    def align(self, srr: str):

        if len(os.listdir(self.config.tools.star.output_folder)) > 0:
            self.log.info(f"STAR output folder is non-empty and assuming "
                          f"star alignment is already finished")
            return

        output_file = f"{self.config.tools.star.output_folder}/{srr}_"

        data = self.meta.data

        opts = [
            "--runThreadN",  # Number of threads for running this program
            str(self.config.no_of_threads),  # If not known keep it as 1
            "--genomeDir",  # Path of the STAR index folder
            self.config.tools.star.index,
            "--outFileNamePrefix",  # Output file prefix
            output_file,
            "--outSAMtype",  # Direct file conversion
            "BAM",  # Directly get BAM file
            "SortedByCoordinate",  # Sort that BAM file with coordinates
            "--readFilesIn"  # Input reads
        ]

        files = self.fasta_for_analysis(srr)
        # Sanity check
        if len(files) not in [1, 2]:
            self.log.error(f"This pipeline currently only supports single "
                           f"end and paired end reads. Your provided files "
                           f"have total of {len(files)} reads")

        analysis_type = "Paired End"
        if len(files) == 1:
            analysis_type = "Single End"

        # Note: Following assumes that forward read in labeled as srr_1
        # while reverse read is labeled as srr_2.
        # TODO: Implement this with configuration file
        files = sorted(files)  # Forward read will come first
        opts.extend(files)

        self.config.tools.star.run(opts, success=f"STAR alignment for {srr} "
                                                 f"is successful.")

        data[srr][self.meta.key_star] = {
            "analysis": analysis_type,
            self.meta.key_input: files,
            "index": self.config.tools.star.index,
            self.meta.key_output: self.config.tools.star.output_folder
        }

        self.meta.append(data)

        # Remove the Log file
        base_file = sys.argv[0]
        log_file = f"{os.path.dirname(base_file)}/Log.out"
        delete_path(log_file)

    def run(self):
        self.index_star()
        for srr in self.meta.get_runs():
            self.align(srr)

        self.log.info("Exiting STAR pipeline")
