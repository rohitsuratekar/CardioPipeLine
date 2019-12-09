#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# Use NCBI SRA tools to manage the SRA files for RNA-seq

import subprocess

from constants.system import SHOW_LOG, FOLDER_NCBI_PUBLIC
from helpers.logging import Log
from helpers.utils import PathUtil


def move_downloaded(sra_id: str, pu: PathUtil, log: Log):
    # Make path for downloaded SRA file
    put_path = f"{pu.sra}/{sra_id}"
    opts = [
        "mkdir",
        "-p",
        put_path
    ]
    if subprocess.run(opts).returncode == 0:
        op2 = [
            "mv",
            pu.ncbi_sra(sra_id),
            pu.sra_file(sra_id)
        ]
        proc = subprocess.run(op2)
        if proc.returncode == 0:
            log.info(f"Downloaded file moved to {put_path}")
        else:
            log.error(f"Error moving file {pu.sra_file(sra_id)} to {put_path}")
    else:
        log.error(f"Error making folder {put_path}")


def fetch_sra(sra_id: str, pu: PathUtil, log: Log):
    opts = [
        pu.prefetch,
        "--verbose",
        sra_id.strip()
    ]
    log.info("Downloading started")
    pros = subprocess.run(opts)
    if pros.returncode == 0:
        log.info(f"{sra_id} is downloaded in {FOLDER_NCBI_PUBLIC}")
        log.info(f"Moving downloaded SRA file to the {pu.sra}")
        move_downloaded(sra_id, pu, log)
    else:
        log.error(f"There was problem downloading {sra_id}")


def get_data(sra_id: str, path_utils: PathUtil, log: Log):
    p = path_utils
    if p.sra_exists(sra_id):
        log.info(f"Raw data is for {sra_id} is already present in "
                 f"{p.sra_file(sra_id)}. Skipping download")
        return
    else:
        log.info(f"Raw file for {sra_id} does not exists in "
                 f"{p.sra_file(sra_id)} and hence fetching from the NCBI "
                 f"server with SRA tools ({p.prefetch})")
        fetch_sra(sra_id, p, log)


def run():
    log = Log(show_log=SHOW_LOG)
    pu = PathUtil()
    get_data("SRX4720625", pu, log)
