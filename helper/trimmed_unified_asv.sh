#!/bin/bash
#SBATCH -p genoa
#SBATCH -t 2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/log/unified_asv_output_%j.out

pipelines="AMPtk VSEARCH DADA2 Tourmaline dadasnake nfcore_ampliseq Natrix DADAist2 QIIME2 QIIME2Deblur MetaWorks APSCALE eDNAFlow LotuS2_otu LotuS2" # 
datasets="mock5 mock12 miseq_v4_sim novaseq_v4_sim gemsim mindset parkinson"  

# datasets=$1 
# echo $datasets

in_dir="/projects/0/gs/danique/pipelines/pipelines"
out_dir="/projects/0/gs/danique/pipelines/int_results"

for d in ${datasets}; do
    for trim_to in 249 199; do
        echo "Processing ${d}"
        echo "Trim to length ${trim_to}"
        
        of="${out_dir}/trim${trim_to}_uni_asv_${d}.fasta"
        if [ -f ${of} ]; then
            rm ${of}
        fi
        touch ${of}

        for p in ${pipelines}; do
            echo "      Processing ${p}"
            if [ $trim_to -eq 249 ] && [[ $p == LotuS2* ]]; then 
                echo "Skipping ${p} for trim ${trim_to} as reads are too short"
                continue
            fi

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
            awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq = seq $0} END {print seq}' $f > tmp.fasta
            f_trim=${f%.fasta}_trim${trim_to}.fasta

            awk -v trim_length="$trim_to" 'BEGIN {FS=">"; OFS=">";} /^>/ {print ">"$2; next} {if(length($0) > trim_length) print substr($0, 1, trim_length); else print $0}' tmp.fasta > ${f_trim}
            f_tmp=${f_trim%.fasta}_prefixed.fasta

            if [ ! -f $f ]; then
                echo "          File ${f} does not exist, skipping"
                continue
            fi

            # add prefix to each sequence header (for tracibility)
            sed "s/^>/>${p}_/" "$f_trim" > "$f_tmp"

            # concatenate to main file
            cat ${f_tmp} >> ${of}

            if [ -f ${f_tmp} ]; then
                rm ${f_tmp}
            fi
            
            d=${d_tmp}

        done 

        # dereplicate (should be 100 match)
        awk '/^>/ {if (seq) print seq; print; seq=""; next} {seq = seq $0} END {print seq}' $of > tmp.fasta
        mv tmp.fasta $of
        uc=${out_dir}/${d}_trim${trim_to}_mapping.uc
        vsearch --derep_fulllength ${of} \
            --output ${of%.fasta}_derep.fasta \
            --uc ${uc} 

    done
done