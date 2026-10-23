#!/bin/bash
#SBATCH -p genoa
#SBATCH -t 2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/log/unified_asv_output_%j.out

# Pools the ASVs of all pipelines for a dataset and dereplicates them with VSEARCH, so the same sequence gets the same ID across pipelines.

# pipelines="AMPtk VSEARCH DADA2 Tourmaline dadasnake nfcore_ampliseq Natrix DADAist2 QIIME2 QIIME2Deblur MetaWorks APSCALE LotuS2 eDNAFlow LotuS2_otu" # NGTax
# datasets="mock5 mock12 miseq_v4_sim novaseq_v4_sim gemsim mindset parkinson"  # mock5 mock12 miseq_v4_sim novaseq_v4_sim

datasets=$1 
echo $datasets

in_dir="/projects/0/gs/danique/pipelines/pipelines"
out_dir="/projects/0/gs/danique/pipelines/int_results"

for d in ${datasets}; do
    echo "Processing ${d}"
    
    of="${out_dir}/uni_asv_${d}.fasta"
    if [ -f ${of} ]; then
        rm ${of}
    fi
    touch ${of}

    for p in ${pipelines}; do
        echo "      Processing ${p}"

        p2="${p%_otu}"
        suffix=""

        if grep -q "_otu" <<< "$p"; then
           suffix="_otu"
        fi  
        
        if [ "$p2" == "Natrix" ] && echo "$d" | grep -q "mock"; then
            d_tmp="${d}"
            d="${d}_o"
        else
            d_tmp="${d}"
        fi

        f=${in_dir}/${p2}/output/${d}/dna-sequences${suffix}_post.fasta
        f_tmp=${f%.fasta}_prefixed.fasta

        if [ ! -f $f ]; then
            echo "          File ${f} does not exist, skipping"
            continue
        fi

        # add prefix to each sequence header (for tracibility)
        sed "s/^>/>${p}_/" "$f" > "$f_tmp"

        # concatenate to main file
        cat ${f_tmp} >> ${of}

        if [ -f ${f_tmp} ]; then
            rm ${f_tmp}
        fi
        
        d=${d_tmp}

    done 

    # dereplicate (should be 100 match)
    uc=${out_dir}/${d}_mapping.uc
    vsearch --derep_fulllength ${of} \
        --output ${of%.fasta}_derep.fasta \
        --uc ${uc} #\
        # --sizeout

    # Add NG-Tax sequences; these sequences seem duplicate (because of lack of merging in pipeline), but are unique (and nonoverlapping with others because of short length)
    f=${in_dir}/NGTax2/output/${d}/dna-sequences_post.fasta
    f_tmp=${f%.fasta}_prefixed.fasta

    if [ ! -f $f ]; then
        echo "          File ${f} does not exist, skipping"
    fi

    # add prefix to each sequence header (for tracibility)
    sed "s/^>/>NGTax2_/" "$f" > "$f_tmp"
    cat "$f_tmp" >> ${of%.fasta}_derep.fasta

    # Add to .uc file
    awk '/^>/ {print substr($0, 2)}' "$f_tmp" | while read -r asv; do
        echo -e "NA\tNA\tNA\tNA\tNA\tNA\tNA\tNA\t${asv}\t*"
    done >> ${uc}

done