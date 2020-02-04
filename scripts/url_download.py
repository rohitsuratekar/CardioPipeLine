#  Copyright (c) 2020, CardioPipeLine
#
#  Author: Rohit Suratekar
#  Institute: Winata Lab, IIMCB Warsaw
#
# Downloads the file from the url file

import urllib.request
from pathlib import Path


def download_from_url():
    Path("downloaded").mkdir(exist_ok=True)
    with open("config/urls.csv") as f:
        for line in f:
            srr, url = line.strip().split(",")
            urllib.request.urlretrieve(url, f"{srr}.sra")
            print(f"{srr} finished ....")
