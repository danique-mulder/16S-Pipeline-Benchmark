#!/bin/bash
#SBATCH -p rome
#SBATCH -t 05:00:00
#SBATCH --ntasks=7
#SBATCH --cpus-per-task=1

module load 2023
module load QIIME2/2023.7

# Demux of all other mocks 
for dir in /gpfs/work2/0/gs/danique/mockrobiota/*/; do 


 fol=$(basename "$dir")
 
 if [ -f "${dir}seq/output/${fol}_mux.qza" ] && [ "$fol" == "mock3" ]; then  
 
 echo "Processing ${fol}"
 
 qiime demux emp-paired \
  --i-seqs "${dir}seq/output/${fol}_mux.qza" \
  --m-barcodes-file "${dir}info/sample-metadata.tsv" \
  --m-barcodes-column BarcodeSequence \
  --p-golay-error-correction \
  --p-rev-comp-barcodes \
  --p-rev-comp-mapping-barcodes \
  --o-per-sample-sequences "${dir}seq/output/${fol}_demux.qza" \
  --o-error-correction-details "${dir}seq/output/${fol}_demux_details.qza" \
  --verbose 
  
 elif [ -f "${dir}seq/output/${fol}_mux.qza" ]; then  
 
 echo "Processing ${fol}"
  
  qiime demux emp-paired \
   --i-seqs "${dir}seq/output/${fol}_mux.qza" \
   --m-barcodes-file "${dir}info/sample-metadata.tsv" \
   --m-barcodes-column BarcodeSequence \
   --p-golay-error-correction \
   --p-rev-comp-mapping-barcodes \
   --o-per-sample-sequences "${dir}seq/output/${fol}_demux.qza" \
   --o-error-correction-details "${dir}seq/output/${fol}_demux_details.qza" \
   --verbose  
   
  else
  
   echo "Skipping ${fol}"
  
   continue
 
 fi
   
 # summarize demultiplex
  qiime demux summarize \
   --i-data "${dir}seq/output/${fol}_demux.qza" \
   --o-visualization "${dir}seq/output/${fol}_demux.qzv"
   
  qiime tools export \
   --input-path "${dir}seq/output/${fol}_demux.qza" \
   --output-path "${dir}fastq_demux"  

done

