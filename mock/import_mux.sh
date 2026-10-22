#!/bin/bash
#SBATCH -p rome
#SBATCH -t 02:00:00
#SBATCH -N 1
#SBATCH --ntasks=7
#SBATCH --cpus-per-task=1

module load 2023
module load QIIME2/2023.7

for dir in /gpfs/work2/0/gs/danique/mockrobiota/*/; do 

 cd $dir
    
 fol=$(basename "$dir")
 
 if [ -f "${dir}seq/barcodes.fastq.gz" ]; then
   qiime tools import \
    --type EMPPairedEndSequences \
    --input-path "${dir}seq" \
    --output-path "${dir}seq/output/${fol}_mux.qza"    
 fi
    
done 


qiime tools import \
    --type EMPPairedEndSequences \
    --input-path "${dir}seq" \
    --output-path "${dir}seq/output/${fol}_mux.qza"
    --imput-format 