#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# Kallisto pipeline


from helpers import ConfigParser, MetaParser, make_path, exists_path
from pipelines._base_pipe import PipeLine


class Kallisto(PipeLine):
    def __init__(self, config: ConfigParser, meta: MetaParser,
                 force_index: bool = False):
        super().__init__()
        self.config = config
        self.meta = meta
        self.sra = meta.sra
        self.force_index = force_index

        out_folder = f"{self.config.names.sra.folder(self.sra)}/kallisto"
        self.config.tools.kallisto.output_folder = out_folder
        make_path(out_folder)
        self.log.info(f"Kallisto output folder set to {out_folder}")

        self.log.info("Kallisto pipeline initialized")

    def index_kallisto(self, index_file: str):
        if exists_path(index_file) and not self.force_index:
            self.log.info(f"Kallisto index is already present at {index_file}")
            self.log.info("If you want to recreate index, "
                          "pass 'force_index=True' as an argument")
            return

        self.log.info("Generating Kallisto index")
        make_path(self.config.tools.kallisto.index)
        opts = [
            "index",  # Index Mode
            "-i",  # Filename (not folder) for indexing
            index_file,
            self.config.names.genome.transcript  # Transcript file
        ]
        self.config.tools.kallisto.run(opts,
                                       success="Kallisto indexing finished "
                                               "successfully.")

    def _single_end(self, srr: str, in_opts: list):
        opts = [x for x in in_opts]
        opts.extend([
            "--single",  # Single end sequencing
            "-l",  # Fragment length
            str(200)  # TODO: Handle this through config file
        ])

        files = self.fasta_for_analysis(srr)

        # Sanity check
        if len(files) != 1:
            self.log.error(
                f"{srr} has total of {len(files)} reads. However, current "
                f"pipeline only supports single end and paried end")

        opts.extend(files)
        self.config.tools.kallisto.run(opts,
                                       success=f"Single end Kallisto quant "
                                               f"analysis on {srr} successful")

    def _paired_end(self, srr: str, in_opts: list):

        opts = [x for x in in_opts]
        opts.extend(self.fasta_for_analysis(srr))
        self.config.tools.kallisto.run(opts,
                                       success=f"Paired end Kallisto quant "
                                               f"analysis on {srr} successful")

    def quant_kallisto(self, index_file: str):

        # Check if quantification is already done
        if exists_path(
                f"{self.config.tools.kallisto.output_folder}/abundance.tsv"):
            self.log.info(
                f"Kallisto quantification is already performed for {self.sra}")
            return

        self.log.info(f"Kallisto quantification for {self.sra} started")
        opts = [
            "quant",  # Quant mode
            "-i",  # Index
            index_file,
            "-o",  # Output folder
            self.config.tools.kallisto.output_folder,
            "-g",  # GTF File (optional)
            self.config.names.genome.gtf,
            "-t",  # Number of threads
            str(self.config.no_of_threads),  # Keep 1 if unknown
            "--bias",  # Sequence based bias correction,
            "--plaintext",  # Plaintext output instead HDF5

        ]

        for srr in self.meta.get_runs():
            data = self.meta.data
            analysis_type = "Paired End"
            if self.meta.is_paired_end(srr):
                self._paired_end(srr, opts)
            else:
                analysis_type = "Single End"
                self._single_end(srr, opts)
            data[srr][self.meta.key_kallisto] = {
                "analysis": analysis_type,
                "input": self.fasta_for_analysis(srr),
                "index": index_file,
                "output": self.config.tools.kallisto.output_folder
            }
            self.meta.append(data)

    def run(self):
        index_file = self.config.tools.kallisto.get_extra("index")
        self.index_kallisto(index_file)
        self.quant_kallisto(index_file)

        self.log.info("Exiting Kallisto pipeline")
