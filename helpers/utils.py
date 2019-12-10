#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# All utils related to this project

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
