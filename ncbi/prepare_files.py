#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# Use NCBI SRA tools to manage the SRA files for RNA-seq

import json
import os
import subprocess
from typing import Union, List

from constants.system import SHOW_LOG, FOLDER_NCBI_PUBLIC, NO_OF_THREADS
from helpers.logging import Log
from helpers.resolver import PathUtil
from helpers.utils import process


def move_downloaded(sra_id: str, pu: PathUtil, log: Log):
    # Make path for downloaded SRA file
    for srr in get_meta_runs(sra_id, pu, log):
        opts = [
            "mv",
            pu.files.ncbi_srr(srr),
            pu.files.srr(sra_id, srr)
        ]
        process(opts,
                log,
                success=f"{srr}.sra moved to {pu.files.srr(sra_id, srr)}",
                error=f"Problem in moving {srr}.sra from "
                      f"{pu.files.ncbi_srr(srr)}"
                )


def fetch_sra(sra_id: str, pu: PathUtil, log: Log):
    opts = [
        pu.tools.prefetch,
        "--verbose",
        sra_id.strip()
    ]
    log.info("Downloading started")

    process(opts, log,
            success=f"{sra_id} runs are downloaded in {FOLDER_NCBI_PUBLIC}",
            error=f"There was problem downloading {sra_id}")

    log.info(f"Moving downloaded SRA runs to the {pu.files.sra_db}")
    move_downloaded(sra_id, pu, log)


def get_meta_runs(sra_id: str, pu: PathUtil, log: Log):
    if pu.sra_meta_exists(sra_id):
        with open(pu.files.sra_meta(sra_id)) as f:
            return json.load(f)["Runs"]
    else:
        pu.make_sra_path(sra_id)
        srr = get_runs(sra_id, pu, log)
        if type(srr) == str:
            srr = [srr]
        all_runs = []
        data = {
            "SRA": sra_id,
        }
        for s in srr:
            current = s.strip().split("\t")
            all_runs.append(current[0])
            data[current[0]] = {
                "Submission": current[1],
                "Status": current[2],
                "Updated": current[3],
                "Published": current[4],
                "Received": current[5],
                "Type": current[6],
                "Center": current[7],
                "Visibility": current[8],
                "Alias": current[9],
                "Experiment": current[10],
                "Sample": current[11],
                "Study": current[12],
                "Loaded": current[13],
                "Spots": current[14],
                "Bases": current[15],
                "Md5sum": current[16],
                "BioSample": current[17],
                "BioProject": current[18],
                "ReplacedBy": current[19]
            }

        data["Runs"] = all_runs
        with open(pu.files.sra_meta(sra_id), 'w', encoding='utf-8') as f:
            json.dump(data, f, indent=4)
            log.info(f"SRA run data is saved in {pu.files.sra_meta(sra_id)}")
        return srr


def prepare_raw_data(sra_id: str, pu: PathUtil, log: Log):
    fetch = False
    for srr in get_meta_runs(sra_id, pu, log):
        if pu.srr_exists(sra_id, srr):
            log.info(f"{srr}.sra run for {sra_id} is already present, "
                     f"skipping download.")
        else:
            fetch = True
    if fetch:
        log.info(f"One or more SRR runs are not present for {sra_id}")
        fetch_sra(sra_id, pu, log)


def get_runs(sra_id: str, pu: PathUtil, log: Log) -> Union[str, List]:
    awk_command = f"{{if ($11==\"{sra_id}\") {{print}} }}"
    opts = [
        "awk",
        "-F",
        " ",
        awk_command,
        pu.sra_accession
    ]
    com = subprocess.run(opts, stdout=subprocess.PIPE)
    com = com.stdout.decode('utf-8').strip().split("\n")
    if len(com) == 0:
        log.error(f"No run data found in the database for {sra_id}. If "
                  f"your SRA id is recently added to the SRA database, "
                  f"please update the SRA_Annotation.tab file and try "
                  f"again. ")
    return com


def convert_to_fastq(sra_id, pu: PathUtil, log: Log, force=False):
    # Check if converted file exists
    fastqs = pu.fastq(sra_id)
    if len(fastqs) == 0 or force:
        log.info("No fastq files found. Converting SRA to fastq")
        for ssr in get_meta_runs(sra_id, pu, log):
            opts = [
                pu.tools.fasterq_dump,  # Tool Path,
                pu.files.srr(sra_id, ssr),  # Path to input SRA file
                "--split-files",  # Splits files according to Reads
                "--progress",  # Show progress
                "--skip-technical",  # Skip technical reads
                "--threads",  # Number of threads
                str(NO_OF_THREADS),  # If you are unsure, use 6 (default)
                "--outdir",  # Output folder
                pu.fastq_dir(sra_id)
            ]
            process(opts, log,
                    success=f"{ssr} is converted into the fastq files",
                    error=f"Error converting {ssr} into fastq")

            all_files = []
            for f in os.listdir(f"{pu.files.sra_db}/{sra_id}"):
                if f.endswith(".fastq") and f.startswith(ssr):
                    all_files.append(f)

            with open(pu.files.sra_meta(sra_id)) as f:
                data = json.load(f)

            data[ssr]["fastq"] = all_files
            with open(pu.files.sra_meta(sra_id), 'w', encoding='utf-8') as f:
                json.dump(data, f, indent=4)

        return pu.fastq(sra_id)
    else:
        log.info(f"Fastq files are already present, skipping conversion")
        log.info(f"If there is problem with converted fastq files, try this "
                 f"function with 'force=True' parameter.")
        return fastqs


if __name__ == "__main()__":
    raise Exception("This script should not be played as a Stand-Alone "
                    "script. Please use it with CardioPipeLine package")


def run():
    log = Log(show_log=SHOW_LOG)
    pu = PathUtil()
    sra_id = "SRX4720626"
    prepare_raw_data(sra_id, pu, log)
    convert_to_fastq(sra_id, pu, log)
