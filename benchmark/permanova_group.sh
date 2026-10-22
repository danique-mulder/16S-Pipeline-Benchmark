#!/bin/bash
#SBATCH -t 2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -o /projects/0/gs/danique/pipelines/log/permanova_group_%A_%a.out
#SBATCH -p rome
#SBATCH --mem=24G
#SBATCH --array=0-4%5

module load 2024
module load R/4.4.2-gfbf-2024a

# Define datasets, levels, and classifiers
datasets=(gemsim mindset parkinson miseq_v4_sim novaseq_v4_sim) 
classifiers=(deftax idtaxa q2 asv) 

# Calculate total number of combinations
n_datasets=${#datasets[@]}
n_classifiers=${#classifiers[@]}
total_combinations=$((n_datasets * n_classifiers))

echo "Total combinations: ${total_combinations}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"

# Get dataset and pipeline for this array task
dataset_idx=$((SLURM_ARRAY_TASK_ID / n_classifiers))
classifier_idx=$((SLURM_ARRAY_TASK_ID % n_classifiers))

# Get names
dataset=${datasets[$dataset_idx]}
classifier=${classifiers[$classifier_idx]}

echo "Processing:"
echo "  Dataset: ${dataset} (index ${dataset_idx})"
echo "  Classifier: ${classifier} (index ${classifier_idx})"

Rscript /projects/0/gs/danique/pipelines/scripts/benchmark/permanova_group.R -d "${dataset}" -c "${classifier}"