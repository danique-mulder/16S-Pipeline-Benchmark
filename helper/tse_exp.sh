#!/bin/bash
#SBATCH -p rome
#SBATCH -t 5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/log/tse_exp_output_%j.out
#SBATCH --mem=200G

# Submits tse_exp.R for each dataset that has a ground truth.

R_LIBS=/projects/0/gs/danique/rpackages:$R_LIBS
export R_LIBS

module load 2024
module load R/4.4.2-gfbf-2024a


datasets="novaseq_v4_sim" # mock5 mock12 miseq_v4_sim 

for d in $datasets; do

    # Run the R script for this specific combination
    Rscript /projects/0/gs/danique/pipelines/scripts/helper/tse_exp.R -d ${d}

done