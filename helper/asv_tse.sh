#!/bin/bash
#SBATCH -t 15:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -o /projects/0/gs/danique/pipelines/log/asv_tse_%A_%a.out
#SBATCH -p rome
#SBATCH --mem-per-cpu=4G
#SBATCH --array=0-27

# SLURM array job that runs unified_asv.sh (to identify identical asv sequences) and asv_tse.R (to combine into a tse) for each dataset × classifier.

module load 2024
module load R/4.4.2-gfbf-2024a

# Define datasets and classifiers
datasets=(mock5 mock12 miseq_v4_sim novaseq_v4_sim gemsim mindset parkinson) 
classifiers=(q2 idtaxa q2_gg idtaxa_gg)  # q2 idtaxa

# Calculate total number of combinations
n_datasets=${#datasets[@]}
n_classifiers=${#classifiers[@]}

# Get dataset and classifier for this array task
dataset_idx=$((SLURM_ARRAY_TASK_ID / n_classifiers))
classifier_idx=$((SLURM_ARRAY_TASK_ID % n_classifiers))

d=${datasets[$dataset_idx]}
c=${classifiers[$classifier_idx]}

echo "Processing dataset: ${d}, classifier: ${c}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"

# Run the unified ASV script once per dataset
if [ $classifier_idx -eq 0 ]; then
    sh /projects/0/gs/danique/pipelines/scripts/helper/unified_asv.sh "${d}"
fi

# Run the R script for this specific combination
Rscript /projects/0/gs/danique/pipelines/scripts/helper/asv_tse.R -d ${d} -c ${c}

