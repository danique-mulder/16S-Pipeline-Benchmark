#!/bin/bash
#SBATCH -t 15:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -o /projects/0/gs/danique/pipelines/log/cor_pip_%A_%a.out
#SBATCH -p rome
#SBATCH --mem=128G
#SBATCH --array=0-6

module load 2024
module load R/4.4.2-gfbf-2024a

# Define datasets, levels, and classifiers
datasets=(mock5 mock12 miseq_v4_sim novaseq_v4_sim gemsim mindset parkinson) # mock5 mock12 miseq_v4_sim novaseq_v4_sim gemsim mindset parkinson
levels=(asv q2 idtaxa q2_gg idtaxa_gg KO cyc asv249 deftax) # asv q2 idtaxa q2_gg idtaxa_gg KO cyc asv249 asv199 asv genus exp KO cyc

# Build combination list
COMBINATIONS=()

for d in "${datasets[@]}"; do
  for lev in "${levels[@]}"; do
  
    if grep -qE "gemsim|mindset|parkinson" <<< $d; then
      e=0
      COMBINATIONS+=("$d $lev $e")
    elif grep -qE "q2$|idtaxa$|KO$|cyc$" <<< "$lev"; then
      e=1
      COMBINATIONS+=("$d $lev $e")
    elif grep -qE "asv$" <<< "$lev"; then
      COMBINATIONS+=("$d $lev 1")
      COMBINATIONS+=("$d $lev 0")
    else 
      e=0
      COMBINATIONS+=("$d $lev $e")
    fi  
  done
done

# Get the total number of combinations
NUM_COMBOS=${#COMBINATIONS[@]}

for i in "${!COMBINATIONS[@]}"; do
  echo "Combination $i: ${COMBINATIONS[$i]}"
done

# Safety check for valid SLURM_ARRAY_TASK_ID
if [[ $SLURM_ARRAY_TASK_ID -ge $NUM_COMBOS ]]; then
  echo "Error: SLURM_ARRAY_TASK_ID ($SLURM_ARRAY_TASK_ID) exceeds available combinations ($NUM_COMBOS)."
  exit 1
fi

# Extract dataset, level, and classifier for this array job
read DATASET LEVEL EXPECTED <<< "${COMBINATIONS[$SLURM_ARRAY_TASK_ID]}"
echo "Running dataset=${DATASET}, level=${LEVEL}, expected=${EXPECTED}"

Rscript /projects/0/gs/danique/pipelines/scripts/benchmark/cor_pip.R -d "${DATASET}" -l "${LEVEL}" -e "${EXPECTED}"