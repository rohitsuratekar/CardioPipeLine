# Converts output from Salmon (which has only transcript IDs) into the genes
#
# Author : Rohit Suratekar
# Date : 21 August 2019


import json
import sys

if len(sys.argv) != 3:
    raise Exception("Missing arguments from the python script. Please pass ["
                    "Path of the annotation.gtf] [path of the salmon's "
                    "quant.sf file]")

GTF_PATH = sys.argv[1]
INPUT_FILE = sys.argv[2]


def extract_attributes(attributes: str):
    attributes = attributes.strip().split(";")
    data = {}
    for a in attributes:
        temp = a.strip().split("\"")
        if len(temp) == 3:
            data[temp[0].strip()] = temp[1].strip()

    return data


def get_annotation():
    with open(GTF_PATH) as f:
        loc = f.tell()
        while f.readline().startswith("#"):
            loc = f.tell()
        f.seek(loc)
        all_data = []
        for line in f:
            data = line.split("\t")
            chromosome = data[0]
            feature = data[2]
            start = data[3]
            end = data[4]
            strand = data[6]
            attributes = extract_attributes(data[8])
            temp = [feature, start, end, chromosome, strand, attributes]
            all_data.append(temp)

        return all_data


def get_data():
    all_data = []
    with open(INPUT_FILE) as f:
        f.readline()  # Skip Header
        for line in f:
            all_data.append(line.strip().split("\t"))

    return all_data


def convert():
    tra = {}
    for d in get_annotation():
        if d[0] == "transcript":
            tra[d[5]["transcript_id"]] = d

    for t in get_data():
        try:
            k = tra[t[0].split(".")[0]]
        except KeyError:
            print(t)


convert()
