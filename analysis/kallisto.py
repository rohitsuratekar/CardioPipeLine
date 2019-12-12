#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# All functions related to Kallisto

from helpers import PathUtil, Log, process


def generate_index(pu: PathUtil, log: Log, force=False):
    # Check if index is already present
    if pu.exists(pu.files.kallisto_index(pu.index.kallisto)) and not force:
        log.info("Kallisto index is already present in the default folder")
        log.info("If you want to regenerate the index, either remove the "
                 "default folder or pass 'force=True'")
        return

    opts = [
        pu.tools.kallisto,  # Path to kallisto binary
        "index",  # Index Mode
        "-i",  # Filename (not folder) for indexing
        pu.files.kallisto_index(pu.index.kallisto),
        pu.transcript_fasta  # Transcript file
    ]

    # Make folder for the index
    pu.make(pu.index.kallisto)

    process(opts, log,
            success="Kallisto index generation is successful",
            error="Something went wrong in generating Kallisto index")


def run():
    pu = PathUtil()
    log = Log(show_log=True)
    generate_index(pu, log)
