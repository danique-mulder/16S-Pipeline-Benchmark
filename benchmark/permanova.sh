#!/bin/bash
#SBATCH -t 1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH -o /projects/0/gs/danique/pipelines/log/permanova_%A_%a.out
#SBATCH -p rome
#SBATCH --mem=64G
#SBATCH --array=0-6

# SLURM array job that runs permanova.R for each dataset × feature level × classifier.

module load 2024
module load R/4.4.2-gfbf-2024a

# Define datasets, levels, and classifiers
datasets=(mock5 mock12 miseq_v4_sim novaseq_v4_sim gemsim mindset parkinson) 
levels=(asv asv249 asv_exp genus exp KO cyc exp_KO exp_cyc ) 
classifiers=(deftax q2 idtaxa idtaxa_gg q2_gg)

# Build combination list
COMBINATIONS=()

for d in "${datasets[@]}"; do
  for lev in "${levels[@]}"; do
    
    if [[ "$lev" == @(asv|asv249|asv199|asv_exp|KO|cyc) ]]; then
      # ASV level or picrust levels: no classifier
      COMBINATIONS+=("$d $lev none")
    
    elif [[ "$lev" == "genus" ]]; then
      # Genus level: both classifiers
      for c in "${classifiers[@]}"; do
        COMBINATIONS+=("$d $lev $c")
      done

    elif [[ "$lev" == "exp" ]]; then
      # Exp level: only for mock* and *_v4_sim datasets, both classifiers
      if [[ "$d" == mock* || "$d" == *_v4_sim ]]; then
        for c in "${classifiers[@]}"; do
          COMBINATIONS+=("$d $lev $c")
        done
      fi
    elif [[ "$lev" == @(exp_KO|exp_cyc) ]] && [[ "$d" == *_v4_sim ]]; then
      COMBINATIONS+=("$d $lev none")
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
read DATASET LEVEL CLASSIFIER <<< "${COMBINATIONS[$SLURM_ARRAY_TASK_ID]}"
echo "Running dataset=${DATASET}, level=${LEVEL}, classifier=${CLASSIFIER}"

# Execute the analysis
if [[ "$CLASSIFIER" == "none" ]]; then
  Rscript /projects/0/gs/danique/pipelines/scripts/benchmark/permanova.R -d "${DATASET}" -l "${LEVEL}"
else
  Rscript /projects/0/gs/danique/pipelines/scripts/benchmark/permanova.R -d "${DATASET}" -c "${CLASSIFIER}" -l "${LEVEL}"
fi

