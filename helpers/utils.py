#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# All utils related to this project

import os
import subprocess

from helpers.logging import Log


def process(options: list, log: Log, *, success: str = None,
            error: str = None) -> bool:
    if subprocess.run(options).returncode == 0:
        if success is not None:
            log.info(success)
        return True
    else:
        if error is not None:
            log.error(error)
        else:
            raise Exception(f"Exception has raised while running following "
                            f"subprocess {options}")
        return False


def get_files_for_mapping(sra_path):
    fastq = []
    for file in os.listdir(sra_path):
        if file.endswith(".fastq"):
            fastq.append(f"{sra_path}/{file}")

    # Check if filtered reads are present
    if os.path.isdir(f"{sra_path}/filtered"):
        filtered = []
        for file in os.listdir(f"{sra_path}/filtered"):
            if file.endswith("filtered.fastq"):
                filtered.append(f"{sra_path}/filtered/{file}")
        if len(filtered) > 0:
            fastq = filtered

    fastq = sorted(fastq)
    return fastq
