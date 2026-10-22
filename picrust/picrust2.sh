#!/bin/bash
#SBATCH -t 12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -o /projects/0/gs/danique/pipelines/log/picrust_%A_%a.out
#SBATCH -p rome
#SBATCH --array=0-118
#SBATCH --mem=224G

CONDA_BASE=$(conda info --base)
CONDA_SH="$CONDA_BASE/etc/profile.d/conda.sh"
source $CONDA_SH 

# dirs
BASE_DIR=/projects/0/gs/danique/pipelines/pipelines
FAILURE_LOG=/projects/0/gs/danique/pipelines/log/picrust_failures.txt

# Specify datasets and pipelines
datasets=(miseq_v4_sim novaseq_v4_sim gemsim mindset parkinson)
pipelines=(AMPtk VSEARCH Tourmaline DADA2 dadasnake NGTax2 nfcore_ampliseq Natrix DADAist2 QIIME2 MetaWorks APSCALE LotuS2 eDNAFlow LotuS2_otu QIIME2Deblur)

# Calculate total number of combinations
n_datasets=${#datasets[@]}
n_pipelines=${#pipelines[@]}
total_combinations=$((n_datasets * n_pipelines))

echo "Total combinations: ${total_combinations}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"

# Get dataset and pipeline for this array task
dataset_idx=$((SLURM_ARRAY_TASK_ID / n_pipelines))
pipeline_idx=$((SLURM_ARRAY_TASK_ID % n_pipelines))

# Get actual names
dataset=${datasets[$dataset_idx]}
pipeline=${pipelines[$pipeline_idx]}

# Function to log failures
log_failure() {
    local reason="$1"
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Task ${SLURM_ARRAY_TASK_ID} FAILED - Dataset: ${dataset}, Pipeline: ${pipeline}, Reason: ${reason}" >> "$FAILURE_LOG"
}

# Safety check
if [ -z "$dataset" ] || [ -z "$pipeline" ]; then
    echo "Invalid array index!"
    echo "  SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
    echo "  dataset_idx: ${dataset_idx}"
    echo "  pipeline_idx: ${pipeline_idx}"
    log_failure "Invalid array index"
    exit 1
fi

echo "Processing:"
echo "  Dataset: ${dataset} (index ${dataset_idx})"
echo "  Pipeline: ${pipeline} (index ${pipeline_idx})"


if [ $pipeline = "LotuS2_otu" ]; then
    seq="${BASE_DIR}/LotuS2/output/${dataset}/dna-sequences_otu_post.fasta"
    feat="${BASE_DIR}/LotuS2/output/${dataset}/feature-table_otu_post.tsv"
    outdir="${BASE_DIR}/LotuS2/output/${dataset}/picrust_otu"
else
    seq="${BASE_DIR}/${pipeline}/output/${dataset}/dna-sequences_post.fasta"
    feat="${BASE_DIR}/${pipeline}/output/${dataset}/feature-table_post.tsv"
    outdir="${BASE_DIR}/${pipeline}/output/${dataset}/picrust"
fi

sed 's/[[:space:]]*$//' "${feat}" > "${feat}.tmp" && mv "${feat}.tmp" "${feat}"

# Check if input files exist
if [ ! -f "$seq" ]; then
    echo "ERROR: FASTA file not found: $seq"
    log_failure "FASTA file not found"
    exit 1
fi

if [ ! -f "$feat" ]; then
    echo "ERROR: Feature table not found: $feat"
    log_failure "Feature table not found"
    exit 1
fi

if [ ! -d "$outdir" ] || [ ! -f "$outdir/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" ] || [ ! -f "$outdir/pathways_out/path_abun_unstrat.tsv.gz" ]; then
    echo "Required files not found, proceeding with Picrust."
else
    echo "Required files already exist, skipping Picrust."
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Task ${SLURM_ARRAY_TASK_ID} EXISTS ALREADY - Dataset: ${dataset}, Pipeline: ${pipeline}" >> "${FAILURE_LOG%.txt}_success.txt"
    exit 0
fi

# Remove output directory if it exists
if [ -d "$outdir" ]; then
    rm -rf "$outdir"
fi


# Run Picrust
conda activate picrust2

picrust2_pipeline.py \
    -s ${seq} \
    -i ${feat} \
    -o ${outdir} \
    -p 16 \
    --max_nsti 2 \
    --verbose \
    --remove_intermediate 

picrust_exit_code=$?
conda deactivate

if [ $picrust_exit_code -ne 0 ]; then
    log_failure "Picrust2 pipeline exited with code $picrust_exit_code"
    exit 1
fi

if [ ! -f "$outdir/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" ] || [ ! -f "$outdir/pathways_out/path_abun_unstrat.tsv.gz" ]; then
    log_failure "Picrust2 pipeline failed - output files not found"
    exit 1
else

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Task ${SLURM_ARRAY_TASK_ID} SUCCEEDED - Dataset: ${dataset}, Pipeline: ${pipeline}" >> "${FAILURE_LOG%.txt}_success.txt"