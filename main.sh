#!/usr/bin/env bash

#bash -e ./preprocessing/requirement.sh

bash -e ./preprocessing/prepare_files.sh
sleep 2m
bash -e ./preprocessing/process_rrna.sh
sleep 2m
# bash ./indexing/star_indexing.sh

bash -e ./alignment/star_alignment.sh
sleep 2m
bash -e ./postprocessing/process_sam_bam.sh
sleep 2m
bash ./postprocessing//clean_files.sh

