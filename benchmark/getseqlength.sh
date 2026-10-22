#!/bin/bash

pipelines="AMPtk VSEARCH Tourmaline DADA2 dadasnake NGTax2 nfcore_ampliseq Natrix DADAist2 QIIME2 MetaWorks APSCALE LotuS2 eDNAFlow LotuS2_otu QIIME2Deblur"
datasets="mock5 mock12 gemsim mindset miseq_v4_sim novaseq_v4_sim parkinson"  

base_dir="/projects/0/gs/danique/pipelines/pipelines"
out_dir="/projects/0/gs/danique/pipelines/int_results"

echo -e "Dataset\tPipeline\tCount\tMin\tMax\tMean\tMedian" > "${out_dir}/seq_length_summary.tsv"
echo -e "Dataset\tPipeline\tID\tLength" > "${out_dir}/seq_length.tsv"

for d in $datasets; do
    for p in $pipelines; do
        echo "Processing ${p} | Dataset ${d}"

        p2=$p
        p=${p%_otu} # without _otu suffix
        suffix=${p2#$p}
        

        if [ ! -f "${base_dir}/${p}/output/${d}/dna-sequences${suffix}_post.fasta" ]; then
            echo "File not found: ${base_dir}/${p}/output/${d}/dna-sequences${suffix}_post.fasta"
            continue
        fi

        d2=$d
        if [[ $d == "mock"* ]] && [[ $p == "Natrix" ]]; then
            d="${d2}_o" 
            echo "Using dataset name ${d} for Natrix (filtering enabled)"
        fi 

        var=$(seqkit fx2tab -nl "${base_dir}/${p}/output/${d}/dna-sequences${suffix}_post.fasta" | awk '{print $2}' | datamash --output-delimiter=$'\t' count 1 min 1 max 1 mean 1 median 1)
        echo -e "${d2}\t${p2}\t${var}" >> "${out_dir}/seq_length_summary.tsv"

        var=$(seqkit fx2tab -nl "${base_dir}/${p}/output/${d}/dna-sequences${suffix}_post.fasta" | awk '{print $1"\t"$2}')
        
        while IFS= read -r line; do
            echo -e "${d2}\t${p2}\t${line}"
        done <<< "$var" >> "${out_dir}/seq_length.tsv"

        d=$d2

    done
done
    


 