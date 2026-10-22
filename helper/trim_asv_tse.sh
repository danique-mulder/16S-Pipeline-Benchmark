#!/bin/bash
#SBATCH -t 15:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -o /projects/0/gs/danique/pipelines/log/asv_tse_%A_%a.out
#SBATCH -p rome
#SBATCH --mem-per-cpu=4G
#SBATCH --array=0-15

module load 2024
module load R/4.4.2-gfbf-2024a

# Define datasets and classifiers
datasets=(miseq_v4_sim novaseq_v4_sim mindset parkinson) # mock5 mock12 gemsim
classifiers=(q2 idtaxa) # q2_gg idtaxa_gg
trim=(249 199) # 

# Calculate dimensions
n_datasets=${#datasets[@]}
n_classifiers=${#classifiers[@]}
n_trim=${#trim[@]}
total_combinations=$((n_datasets * n_classifiers * n_trim))

echo "Total combinations: $total_combinations (datasets: $n_datasets, classifiers: $n_classifiers, trim: $n_trim)"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"

# dataset × classifier × trim
trim_idx=$(( SLURM_ARRAY_TASK_ID % n_trim ))
classifier_idx=$(( (SLURM_ARRAY_TASK_ID / n_trim) % n_classifiers ))
dataset_idx=$(( SLURM_ARRAY_TASK_ID / (n_trim * n_classifiers) ))

d=${datasets[$dataset_idx]}
c=${classifiers[$classifier_idx]}
t=${trim[$trim_idx]}

echo "Processing: dataset=${d}, classifier=${c}, trim=${t}"
echo "Indices: dataset=$dataset_idx, classifier=$classifier_idx, trim=$trim_idx"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"

# Run the unified ASV script once per dataset
if [ $trim_idx -eq 0 ] && [ $classifier_idx -eq 0 ]; then
    echo "Running ASV preprocessing for dataset: ${d}"
    sh /projects/0/gs/danique/pipelines/scripts/helper/trimmed_unified_asv.sh "${d}"
fi

# Run the R script for this specific combination
Rscript /projects/0/gs/danique/pipelines/scripts/helper/asv_tse.R -d ${d} -c ${c} -t ${t}