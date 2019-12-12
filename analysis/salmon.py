#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# All Salmon related analysis

from constants.system import NO_OF_THREADS
from helpers import PathUtil, Log, process, get_files_for_mapping


def generate_decoy(pu: PathUtil, log: Log):
    # Check if file exists already
    if pu.exists(pu.files.salmon_decoy_keys(pu.index.salmon)):
        log.info("Salmon decoy index files are already available, skipping "
                 "their generation.")
        return
    log.info("Salmon decoy file generation started.")
    opts = [
        pu.tools.salmon_decoy,  # Path to decoy script
        "-j",  # Number of threads
        str(NO_OF_THREADS),  # If unknown, keep 1 or remove option
        "-a",  # GTF genome annotation
        pu.gtf_annotation,
        "-g",  # Genome in fasts
        pu.genome_fasta,
        "-t",  # Transcriptom from cDNA in fasta
        pu.transcript_fasta,
        "-o",  # Output path
        pu.index.salmon,  # Keep in index folder itself
        "-b",  # Path to bedtools binary
        pu.tools.bedtools,
        "-m",  # Path to mashmap binary
        pu.tools.mashmap
    ]

    process(opts, log,
            success=f"Salmon Decoy file generated successfully in "
                    f"{pu.index.salmon}",
            error="There was problem in generating Salmon decoy file")


def generate_index(pu: PathUtil, log: Log, force=False):
    # Quickly check if index is present by checking info.json

    if pu.exists(pu.files.salmon_index_info(pu.index.salmon)) and not force:
        log.info("Salmon Index's info.json file is present in the default "
                 "index folder and hence assuming indexing is done")
        log.info("If you want to reindex, use 'force=True' option")
        return

    # Note : Use decoy index generated by above function instead of original
    # transcript for indexing
    # Relevant information is on following tutorial
    # https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
    opts = [
        pu.tools.salmon,  # Salmon binary
        "index",  # Index mode
        "-t",  # Decoy transcript path (Read above note)
        pu.files.salmon_decoy_index(pu.index.salmon),
        "-i",  # Index folder
        pu.index.salmon,
        "--decoys",  # Point to Decoy Text file
        pu.files.salmon_decoy_keys(pu.index.salmon),
        "-k",  # k-mer length.
        str(31),  # 31 is default and recommended by the authors
        "-p",  # Number of threads
        str(NO_OF_THREADS)  # If unknown, don't use it
    ]

    # Make the index folder
    pu.make(pu.index.salmon)

    log.info("Salmon indexing started")

    process(opts, log,
            success="Salmon indexing is successfully finished",
            error="Something went wrong in generating Salmon index")


def salmon_quant(sra_id: str, pu: PathUtil, log: Log):
    opts = [
        pu.tools.salmon,  # Salmon binary
        "quant",  # Quant mode
        "-i",  # Index path
        pu.index.salmon,
        "-l",  # Library type
        "A",  # Automatic
        "--validateMapping",
        # enables selective alignment of the sequencing reads
        "--gcBias",  # Corrects for  fragment-level GC biases. This option
        # is also recommended by DESeq2 analysis
        "-p",  # Number of threads
        str(NO_OF_THREADS),  # Remove this option if not sure
        "-g",  # Additional GTF file to get quants.genes.sf
        pu.gtf_annotation,
        "-o",  # Output path
        pu.files.salmon_output_folder(sra_id, make=True)  # If not present,
        # create the output folder path with 'make=True' argument
    ]

    # Add files to the options
    fastq_files = get_files_for_mapping(pu.files.sra_folder(sra_id))

    if not len(fastq_files) == 2:
        log.error("Single end sequencing pipeline is not set yet.")

    for i, f in enumerate(fastq_files):
        opts.append(str(-(i + 1)))
        opts.append(f)

    log.info("Starting Salmon quantification")
    process(opts, log,
            success=f"Salmon quantification for {sra_id} is successful",
            error=f"Something went wrong with {sra_id} Salmon quantification")


def run():
    pu = PathUtil()
    log = Log(show_log=True)
    sra_id = "SRX4720626"
    generate_decoy(pu, log)
    generate_index(pu, log)
    salmon_quant(sra_id, pu, log)