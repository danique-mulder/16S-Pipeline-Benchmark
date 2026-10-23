#!/bin/bash

# Renames the GEMSIM FASTQ files to a standard naming scheme, writes a manifest that links new and old names, and submits QC.

DATA_DIR="/gpfs/work2/0/gs/raw_data/gemsim/16s"
NEW_DATA_DIR="/gpfs/work2/0/gs/danique/pipelines/data/gemsim"

mkdir -p $NEW_DATA_DIR
mkdir -p "${NEW_DATA_DIR}/seq"

# Write the header row to the manifest file
echo -e "new_name\told_name\tid\tbatch\tposition\tbarcode\tlane\tread\tset\tflowid" > "${NEW_DATA_DIR}/manifest.txt"

for file in ${DATA_DIR}/*.fastq.gz; do

    old_name=$(basename $file)

    # Split the filename by underscores
    IFS='_' read -r -a parts <<< "$old_name"

    # Remove hyphens from the first part
    sample=$(echo "${parts[0]}" | tr -d '-')
    file_extension="${parts[3]}_${parts[4]}_${parts[5]}"
    batch=$(echo "${parts[1]}" | awk -F'-' '{print $1}')
    position=$(echo "${parts[1]}" | awk -F'-' '{print $2}')
    barcode="${parts[2]}"
    lane="${parts[3]}"
    read="${parts[4]}"
    set="${parts[5]}"
    flowid="${parts[6]}"

    # Construct the new filename
    new_name="${sample}_${file_extension}.fastq.gz"
    echo $new_name
    # add to manifest
    echo -e "${new_name}\t${old_name}\t${sample}\t${batch}\t${position}\t${barcode}\t${lane}\t${read}\t${set}\t${flowid}" >> "${NEW_DATA_DIR}/manifest.txt"

    ln -sf $file "${NEW_DATA_DIR}/seq/${new_name}"
done

head "${NEW_DATA_DIR}/manifest.txt"

QC=/projects/0/gs/danique/pipelines/scripts/qc/qc.sh
sbatch -t 01:00:00 $QC -d "gemsim" -i "${NEW_DATA_DIR}/seq"



