#!/bin/bash
#SBATCH -t 5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH -o /projects/0/gs/danique/pipelines/log/picrust_sim%j.out
#SBATCH -p rome

CONDA_BASE=$(conda info --base)
CONDA_SH="$CONDA_BASE/etc/profile.d/conda.sh"
source $CONDA_SH 

dataset=sim
seq=/projects/0/gs/danique/pipelines/simulations/input/sequences_v4.fasta
feat=/projects/0/gs/danique/pipelines/simulations/input/implanted_counts.tsv
outdir=/projects/0/gs/danique/pipelines/simulations/picrust 

if [ ! -d "$outdir" ] || [ ! -f "$outdir/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz" ] || [ ! -f "$outdir/pathways_out/path_abun_unstrat.tsv.gz" ]; then
    echo "Required files not found, proceeding with Picrust."
else
    echo "Required files already exist, skipping Picrust."
    exit 0
fi

if [ -d "$outdir" ]; then
    rm -rf "$outdir"
fi

# Run Picrust
conda activate picrust2

echo "Start picrust"

picrust2_pipeline.py \
    -s ${seq} \
    -i ${feat} \
    -o ${outdir} \
    -p 16 \
    --max_nsti 2 \
    --verbose \
    --remove_intermediate

echo "Picrust completed successfully"

conda deactivate