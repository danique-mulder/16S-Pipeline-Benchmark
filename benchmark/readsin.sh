#!/bin/bash
#SBATCH -p rome
#SBATCH -t 6:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/gpfs/work2/0/gs/danique/test/readsin_%j.out

datasets="mock5 mock12 gemsim mindset parkinson miseq_v4_sim novaseq_v4_sim" 

for d in $datasets; do
    echo "Processing dataset: $d"

    if echo "$d" | grep -q mock; then
        datapath="/projects/0/gs/danique/mockrobiota/$d/fastq_demux"
    elif echo "$d" | grep -q gemsim; then
        datapath="/projects/0/gs/danique/pipelines/data/gemsim/seq"
    elif echo "$d" | grep -q mindset; then
        datapath="/projects/0/gs/danique/pipelines/data/mindset/seq"
    elif echo "$d" | grep -q parkinson; then
        datapath="/projects/0/gs/danique/pipelines/data/parkinson/seq"
    elif echo "$d" | grep -q miseq; then
        datapath="/projects/0/gs/danique/pipelines/simulations/input/fastq/miseq"
    elif echo "$d" | grep -q novaseq; then
        datapath="/projects/0/gs/danique/pipelines/simulations/input/fastq/novaseq"
    fi

    sample_names=$(head /projects/0/gs/danique/pipelines/pipelines/DADA2/output/$d/feature-table_post.tsv -n 1 | cut -f 2- | tr '\t' '\n')
    
    output_file="/projects/0/gs/danique/pipelines/int_results/routput/${d}.tsv"
    echo -e "ID\treads_in" > "$output_file"

    for sample in $sample_names; do

        echo "Processing sample: $sample"
        file=$(ls ${datapath} | grep "${sample}" | grep "R1" | head -n 1)
        file="${datapath}/${file}"

        if [ -n "$file" ]; then
            count=$(seqtk seq -A "$file" | grep -c "^>")
            echo -e "${sample}\t${count}" >> "$output_file"
        else
            echo "File for sample $sample not found in $datapath"
        fi
    done
done
