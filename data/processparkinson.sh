#!/bin/bash

# Renames the Parkinson FASTQ files to a standard naming scheme, writes a manifest that links new and old names, and submits QC.


DATA_DIR="/gpfs/work2/0/gs/danique/raw_data/PRJNA601994/16s"
NEW_DATA_DIR="/gpfs/work2/0/gs/danique/pipelines/data/parkinson"

mkdir -p "${NEW_DATA_DIR}/seq"

# Write the header row to the manifest file
echo -e "new_name\told_name\tid\tread" > "${NEW_DATA_DIR}/manifest.txt"

for file in ${DATA_DIR}/*.fastq.gz; do

    old_name=$(basename $file)

    # Split the filename by underscores
    IFS='_' read -r -a parts <<< "$old_name"

    # Remove hyphens from the first part
    sample=$(echo "${parts[0]}" | tr -d '-')
    read="${parts[1]%.fastq.gz}"

    new_name="${sample}_L001_R${read}_001.fastq.gz"

    # Construct the new filename
    echo $new_name
 
    # add to manifest
    echo -e "${new_name}\t${old_name}\t${sample}\t${read}" >> "${NEW_DATA_DIR}/manifest.txt"

    ln -sf $file "${NEW_DATA_DIR}/seq/${new_name}"
done

head "${NEW_DATA_DIR}/manifest.txt"

# QC
QC=/projects/0/gs/danique/pipelines/scripts/qc/qc.sh
sbatch -t 01:00:00 $QC -d "parkinson" -i "${NEW_DATA_DIR}/seq"


# Create combined metadata file 
Rscript - <<EOF
.libPaths(c("/projects/0/gs/danique/rpackages", .libPaths()))
library(dplyr)
library(readr)

meta1 <- "/projects/0/gs/danique/raw_data/PRJNA601994/sampleid.tsv"
meta2 <- "/projects/0/gs/danique/raw_data/PRJNA601994/PRJNA601994_key.tsv"
meta3 <- "/projects/0/gs/danique/raw_data/PRJNA601994/PRJNA601994_sample_attributes_HP.tsv"
meta_out <- "/projects/0/gs/danique/raw_data/metadata.tsv"

m1 <- read_tsv(meta1) 
m2 <- read_tsv(meta2)
m3 <- read_tsv(meta3)
colnames(m3)[1] <- "sample_accession"

meta_all <- m1 %>%
  dplyr::left_join(m2, by = "run_accession") %>%
  dplyr::left_join(m3, by = "sample_accession")

colnames(meta_all)

write.table(meta_all, meta_out, row.names = FALSE, quote = FALSE, sep = "\t")
EOF

cp /projects/0/gs/danique/raw_data/metadata.tsv /projects/0/gs/danique/pipelines/int_results/metadata/parkinson.tsv
