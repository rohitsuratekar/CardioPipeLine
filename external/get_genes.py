# Converts output from Salmon (which has only transcript IDs) into the genes
#
# Author : Rohit Suratekar
# Date : 21 August 2019
#
# Ensemble annotation file can be downloaded from the BioMart with following
# columns
#
# Gene stable ID,	Transcript stable ID,	Chromosome/scaffold
# name,	Gene start (bp),	Gene end (bp),	Strand,	Transcription start site (
# TSS),	Transcript name,	Gene name,	Source (gene),	Source (
# transcript),Transcript type,	Transcript length (including UTRs and CDS)	,
# ZFIN ID

import sys

import pandas as pd


def get_data_frame(input_file: str):
    """
    Get data from the file and converts into the pandas DataFrame
    :param input_file: Full path of the file
    :return: Pandas DataFrame
    """
    with open(input_file) as f:
        return pd.read_csv(f, delimiter="\t", low_memory=False)


def convert(input_file: str, annotation_file: str, output_file: str):
    """
    Main function which takes the Ensemble annotation file (derived from
    BioMart) and joins it with Salmon output file
    :param input_file: Salmon output file
    :param annotation_file: Ensemble annotation file
    :param output_file: Path of the out file
    """

    tra_id = "MainTranscript"
    ann_tra_id = "Transcript stable ID"
    ann_gene_id = "Gene stable ID"
    ann_gene_name = "Gene name"
    data_name = "Name"
    data = get_data_frame(input_file)
    ann = get_data_frame(annotation_file)
    data[tra_id] = data[data_name].str.split(".").str[0]
    ann = ann[[ann_gene_id, ann_tra_id, ann_gene_name, "ZFIN ID",
               "Source (transcript)"]]
    data = data.merge(ann, left_on=tra_id, right_on=ann_tra_id, how="inner")
    data = data[~data.duplicated()]
    del data[tra_id]
    del data[ann_tra_id]

    with open(output_file, "w") as f:
        print(data.round(4).to_csv(sep="\t", index=False), file=f)


if __name__ == "__main__":

    if len(sys.argv) != 4:
        raise Exception(
            "Missing arguments from the python script. Please pass ["
            "Path of the transcript_to_gene.tsv] [path of the "
            "salmon's quant.sf file] [path of output file]")

    ANN_FILE = sys.argv[1]
    INPUT_FILE = sys.argv[2]
    OUTPUT_FILE = sys.argv[3]

    convert(INPUT_FILE, ANN_FILE, OUTPUT_FILE)
