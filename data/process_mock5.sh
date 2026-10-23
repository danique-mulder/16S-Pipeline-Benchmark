#!/bin/bash

# Imports data into QIIME2, demultiplexed and addds primers to the reads

# # Processing on mock5

module load 2023
module load QIIME2/2023.7

# EMP barcodes, ran from /gpfs/work2/0/gs/danique/mockrobiota/
cd /gpfs/work2/0/gs/danique/mockrobiota/

# Import multiplexed sequences for mock5
if [ -f "mock5/seq/barcodes.fastq.gz" ]; then
qiime tools import \
   --type EMPPairedEndSequences \
   --input-path "mock5/seq" \
   --output-path "mock5/seq/output/mock5_mux.qza"    
else 
   echo "Barcodes file not found in mock5/seq. Skipping import."
   exit 1
fi

# demux
qiime demux emp-paired \
   --i-seqs "mock5/seq/output/mock5_mux.qza" \
   --m-barcodes-file "mock5/info/sample-metadata.tsv" \
   --m-barcodes-column BarcodeSequence \
   --p-golay-error-correction \
   --p-rev-comp-mapping-barcodes \
   --o-per-sample-sequences "mock5/seq/output/mock5_demux.qza" \
   --o-error-correction-details "mock5/seq/output/mock5_demux_details.qza" \
   --verbose 

# summarize demux
qiime demux summarize \
   --i-data "mock5/seq/output/mock5_demux.qza" \
   --o-visualization "mock5/seq/output/mock5_demux.qzv" 
    
# Export as fastq files
qiime tools export \
   --input-path "mock5/seq/output/mock5_demux.qza" \
   --output-path "mock5/fastq_demux"

# Processing on mock5
BASE_DIR="/gpfs/work2/0/gs/danique/mockrobiota"
DATA_DIR="${BASE_DIR}/mock5"

# Rename files to <sample_name>_L001_R[12]_001.fastq.gz
for file in ${DATA_DIR}/fastq_demux/*.fastq.gz; do
    new_name=$(basename "$file" | sed 's/HMPMockV1.[12].//; s/_[1234]//')
    echo $new_name
    mv "$file" "${DATA_DIR}/fastq_demux/$new_name"
done

# Add primers
SCRIPT_DIR="/gpfs/work2/0/gs/danique/pipelines/scripts/helper"
sh ${SCRIPT_DIR}/addprimers.sh ${DATA_DIR}/fastq_demux

# Added primers are: 
#forward="GTGCCAGCAGCCGCGGTAA"
#reverse="GGACTACAAGGGTATCTAAT"

# remultiplex sequences
sh ${SCRIPT_DIR}/remux.sh -i ${DATA_DIR}/fastq_demux -p ${DATA_DIR}/fastq_demux/primer -f GTGCCAGCAGCCGCGGTAA -r GGACTACAAGGGTATCTAAT




