#!/bin/bash
#SBATCH -p fat_rome
#SBATCH -t 10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/log/tse_filter_output_%j.out
#SBATCH --mem=500G

# Submits tse_filter.R

module load 2024
module load R/4.4.2-gfbf-2024a
Rscript /projects/0/gs/danique/pipelines/scripts/benchmark/tse_filter.R