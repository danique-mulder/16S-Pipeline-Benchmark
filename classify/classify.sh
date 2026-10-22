#!/bin/bash

# Specify pipelines and datasets to run classification on
pipelines="NGTax2 VSEARCH AMPtk DADA2 Tourmaline dadasnake nfcore_ampliseq Natrix QIIME2 QIIME2Deblur MetaWorks APSCALE LotuS2 DADAist2 eDNAFlow" 
datasets="mock5 mock12 gemsim mindset parkinson novaseq_v4_sim miseq_v4_sim" 

# Classification sciprts
q2_nb="/gpfs/work2/0/gs/danique/pipelines/scripts/classify/q2_naivebayes.sh" # script for q2 naive bayes
idtaxa="/gpfs/work2/0/gs/danique/pipelines/scripts/classify/IDTAXA.sh" # script for IDTAXA

LOG_FILE="/gpfs/work2/0/gs/danique/pipelines/log/log_classify.tsv"
echo -e "classifier\tjobid\tpipeline\tdataset" > "${LOG_FILE}"

for dataset in $datasets; do
    for pipeline in $pipelines; do
        for suffix in "" "_otu"; do   # "" "otu"
        
            if [[ $suffix == "_otu" ]] && [[ $pipeline != "LotuS2" ]]; then
                    continue
            fi

        seq="dna-sequences${suffix}_post.fasta"
        echo $pipeline
        echo $seq

        echo "Running Q2 naive bayes for $pipeline on $dataset with sequence file $seq"
        command="sbatch -p genoa $q2_nb -d $dataset -p $pipeline -s $seq"
        jobid=$($command | awk '{print $4}')
        echo -e "q2\t${jobid}\t${pipeline}\t${dataset}" >> "${LOG_FILE}" 
        echo 

        echo "Running IDTAXA for $pipeline on $dataset with sequence file $seq"
        command="sbatch -p genoa $idtaxa -d $dataset -p $pipeline -s $seq"
        jobid=$($command | awk '{print $4}')
        echo -e "idtaxa\t${jobid}\t${pipeline}\t${dataset}" >> "${LOG_FILE}" 

        # If pipeline == NGTax, run classification on the concatenated forward and reverse sequences too
        if [ $pipeline == "NGTax2" ]; then
        
            seq="dna-sequences_concat.fasta"

            echo "Running Q2 naive bayes for $pipeline on $dataset with sequence file $seq"
            command="sbatch -p genoa $q2_nb -d $dataset -p $pipeline -s $seq"
            jobid=$($command | awk '{print $4}')
            echo -e "q2\t${jobid}\t${pipeline}\t${dataset}" >> "${LOG_FILE}" 
            echo 

            echo "Running IDTAXA for $pipeline on $dataset with sequence file $seq"
            command="sbatch -p genoa $idtaxa -d $dataset -p $pipeline -s $seq"
            jobid=$($command | awk '{print $4}')
            echo -e "idtaxa\t${jobid}\t${pipeline}\t${dataset}" >> "${LOG_FILE}" 
        fi

        done
    done
done

