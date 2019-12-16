#   Copyright (c)  2020, CardioPipeLine
#   Author: Rohit Suratekar
#   Organization: IIMCB
#
# Util functions

import os
import shutil


def exists_path(path: str) -> bool:
    """
    Checks if given file or folder exists
    :param path: Path of the file or folder
    :return: True, if it exists
    """
    return os.path.exists(path)


def make_path(path: str):
    """
    Simple function to make folder structure
    :param path: Path which you need to make
    """
    if os.path.exists(path):
        # If it is already exists, just return
        return
    # make the folder structure
    os.makedirs(path)


def delete_path(path: str):
    """
    Simple context aware tool to remove file or folder
    :param path: Path of the file or folder
    """
    if not os.path.exists(path):
        # If path doesn't exist, just return
        return
    if os.path.isfile(path):
        # If it is a file
        os.remove(path)
    else:
        # If it is a folder
        shutil.rmtree(path)


def move_path(origin: str, destination: str):
    """
    Moves files or folder from origin to destination
    :param origin: Path of the origin file or folder
    :param destination: Path of the destination file or folder
    """
    shutil.move(origin, destination)
