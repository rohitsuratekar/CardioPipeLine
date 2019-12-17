#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# Functions related to StringTie

from helpers import ConfigParser, MetaParser, make_path, exists_path
from pipelines._base_pipe import PipeLine


class StringTie(PipeLine):
    def __init__(self, config: ConfigParser, meta: MetaParser):
        super().__init__()
        self.config = config
        self.meta = meta
        self.sra = meta.sra

        out_folder = f"{self.config.names.sra.folder(self.sra)}/stringtie"
        make_path(out_folder)
        self.config.tools.stringtie.output_folder = out_folder
        self.log.info(f"StringTie output folder set as {out_folder}")

        self.log.info("StringTie pipeline initialized")

    def get_input_file(self, srr):
        return self.meta.data[srr][self.meta.key_star]["output"]

    def quantify(self, srr: str):
        output_gtf = f"{self.config.tools.stringtie.output_folder}/{srr}_assembled.gtf"
        output_file = f"{self.config.tools.stringtie.output_folder}/{srr}_gene_expression.tsv"
        bam_file = f"{self.get_input_file(srr)}/{srr}_Aligned.sortedByCoord.out.bam"

        if exists_path(output_file):
            self.log.info(f"StringTie analysis is already performed on {srr}")
            return

        opts = [
            "-p",  # Number of threads
            str(self.config.no_of_threads),
            "-G",  # Annotation
            self.config.names.genome.gtf,
            "-e",
            # only estimate the abundance of given reference transcripts (
            # requires -G). You might want to ommit this if working with
            # novel transcripts
            "-B",  # enable output of Ballgown table files
            "-o",  # Output GTF file name
            output_gtf,
            "-A",  # Gene abundances file
            output_file,
            bam_file  # Input file
        ]
        self.config.tools.stringtie.run(opts,
                                        success=f"StringTie quantification "
                                                f"for {srr} successful ")

        # Write metadata

        data = self.meta.data
        data[srr][self.meta.key_string_tie] = {
            "input": bam_file,
            "output": output_file
        }
        self.meta.append(data)

    def run(self):
        for srr in self.meta.get_runs():
            self.quantify(srr)

        self.log.info("Exiting StringTie pipeline")
