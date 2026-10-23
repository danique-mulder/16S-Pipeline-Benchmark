#!/bin/bash
#SBATCH -p rome
#SBATCH -t 1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/log/runrmd_output_%j.out
#SBATCH --mem=64G

# Renders pipelines.Rmd to HTML.

R_LIBS=/projects/0/gs/danique/rpackages:$R_LIBS
export R_LIBS

module load 2024
module load R/4.4.2-gfbf-2024a

Rscript -e "rmarkdown::render('/projects/0/gs/danique/pipelines/scripts/benchmark/pipelines.Rmd', 'html_document', output_dir = '/projects/0/gs/danique/pipelines/int_results/routput')"

