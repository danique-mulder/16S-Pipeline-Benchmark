#!/bin/bash

# processes the mock12 dataset and adds primers to reads 

# Since this mock community only consists of one sample, there is no need for demultiplexing the reads.
BASE_DIR="/gpfs/work2/0/gs/danique/mockrobiota"
DATA_DIR="${BASE_DIR}/mock12"

zcat "${DATA_DIR}/seq/forward.fastq.gz" | \
 awk 'NR % 4 == 1 {if ($1 ~ /\.[12]$/) sub(/\.[12]$/, "", $1);} NR % 4 == 3 {if ($1 ~ /\.[12]$/) sub(/\.[12]$/, "", $1);} {print}' | \
 gzip > "${DATA_DIR}/fastq_demux/mock12_L001_R1_001.fastq.gz"

zcat "${DATA_DIR}/seq/reverse.fastq.gz" | \
 awk 'NR % 4 == 1 {if ($1 ~ /\.[12]$/) sub(/\.[12]$/, "", $1);} NR % 4 == 3 {if ($1 ~ /\.[12]$/) sub(/\.[12]$/, "", $1);} {print}' | \
 gzip > "${DATA_DIR}/fastq_demux/mock12_L001_R2_001.fastq.gz"

# Add primers
SCRIPT_DIR="/gpfs/work2/0/gs/danique/pipelines/scripts/helper"
sh ${SCRIPT_DIR}/addprimers.sh ${DATA_DIR}/fastq_demux



