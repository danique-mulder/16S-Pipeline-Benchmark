#!/bin/bash

DATASET=$1
STAT_FILE="/gpfs/work2/0/gs/danique/pipelines/log/generalStats_${DATASET}.tsv"
LOG_FILE="/gpfs/work2/0/gs/danique/pipelines/log/checkformat_${DATASET}.tsv"

# Initialize the stats file
echo -e "pipeline\tnumSeq\tnumFeat\tnumGenus\tnumGenus_idTaxa\tnumgenus_q2\tnumGenus_idTaxa_gg\tnumgenus_q2_gg" > "$STAT_FILE"
touch "$LOG_FILE"

# Define pipelines and their types
declare -A PIPELINES=(
    ["ASV"]="DADA2 LotuS2 QIIME2 QIIME2Deblur VSEARCH AMPtk DADAist2 dadasnake eDNAFlow Natrix nfcore_ampliseq Tourmaline APSCALE MetaWorks NGTax2"
    ["OTU"]="LotuS2"
)

# Function to process a pipeline
process_pipeline() {
    local pipeline=$1
    local type=$2
    local suffix=$3

    echo "===== Pipeline: $pipeline ====="
    echo "===== Pipeline: $pipeline =====" >> "$LOG_FILE"

    echo "dna-sequences.fasta:" >> "$LOG_FILE"
    
    if [[ -f "/gpfs/work2/0/gs/danique/pipelines/pipelines/${pipeline}/output/${DATASET}/dna-sequences${suffix}_post.fasta" ]]; then
        head -n 2 "/gpfs/work2/0/gs/danique/pipelines/pipelines/${pipeline}/output/${DATASET}/dna-sequences${suffix}_post.fasta" >> "$LOG_FILE"
    else
        echo "/gpfs/work2/0/gs/danique/pipelines/pipelines/${pipeline}/output/${DATASET}/dna-sequences${suffix}_post.fasta does not exist" >> "$LOG_FILE"
    fi

    echo >> "$LOG_FILE"
    echo "feature-table.tsv:" >> "$LOG_FILE"

    if [[ -f "/gpfs/work2/0/gs/danique/pipelines/pipelines/${pipeline}/output/${DATASET}/dna-sequences${suffix}_post.fasta" ]]; then
        head -n 2 "/gpfs/work2/0/gs/danique/pipelines/pipelines/${pipeline}/output/${DATASET}/feature-table${suffix}_post.tsv" >> "$LOG_FILE"
    else
        echo "/gpfs/work2/0/gs/danique/pipelines/pipelines/${pipeline}/output/${DATASET}/feature-table${suffix}_post.tsv does not exist" >> "$LOG_FILE"
    fi

    echo >> "$LOG_FILE"
    echo "taxonomy.tsv:" >> "$LOG_FILE"

    if [[ -f "/gpfs/work2/0/gs/danique/pipelines/pipelines/${pipeline}/output/${DATASET}/taxonomy${suffix}_post.tsv" ]]; then
        head -n 2 "/gpfs/work2/0/gs/danique/pipelines/pipelines/${pipeline}/output/${DATASET}/taxonomy${suffix}_post.tsv" >> "$LOG_FILE"
    else
        echo "/gpfs/work2/0/gs/danique/pipelines/pipelines/${pipeline}/output/${DATASET}/taxonomy${suffix}_post.tsv does not exist" >> "$LOG_FILE"
    fi
    echo >> "$LOG_FILE"

    # Calculate numasv and numgenus
    local numseq=$(awk '/^>/' "/gpfs/work2/0/gs/danique/pipelines/pipelines/${pipeline}/output/${DATASET}/dna-sequences${suffix}_post.fasta" | wc -l)
    local numfeat=$(sed '1d' "/gpfs/work2/0/gs/danique/pipelines/pipelines/${pipeline}/output/${DATASET}/feature-table${suffix}_post.tsv" | wc -l)
    local numgenus=$(sed '1d' "/gpfs/work2/0/gs/danique/pipelines/pipelines/${pipeline}/output/${DATASET}/taxonomy${suffix}_post.tsv" | awk -F'\t' '{print $7}' | sort | uniq | wc -l)
    local numgenus_idtaxa=$(sed '1d' "/gpfs/work2/0/gs/danique/pipelines/pipelines/${pipeline}/output/${DATASET}/idtaxa${suffix}_V4.tsv" | awk -F'\t' '{print $7}' | sort | uniq | wc -l)
    local numgenus_q2=$(sed '1d' "/gpfs/work2/0/gs/danique/pipelines/pipelines/${pipeline}/output/${DATASET}/q2_nb_silva${suffix}.tsv" | awk -F'\t' '{print $7}' | sort | uniq | wc -l)
    local numgenus_idtaxa_gg=$(sed '1d' "/gpfs/work2/0/gs/danique/pipelines/pipelines/${pipeline}/output/${DATASET}/idtaxa${suffix}_gg_V4.tsv" | awk -F'\t' '{print $7}' | sort | uniq | wc -l)
    local numgenus_q2_gg=$(sed '1d' "/gpfs/work2/0/gs/danique/pipelines/pipelines/${pipeline}/output/${DATASET}/q2_nb_gg${suffix}.tsv" | awk -F'\t' '{print $7}' | sort | uniq | wc -l)


    # Append to stats file
    echo -e "${pipeline}${suffix}\t${numseq}\t${numfeat}\t${numgenus}\t${numgenus_idtaxa}\t${numgenus_q2}\t${numgenus_idtaxa_gg}\t${numgenus_q2_gg}" >> "$STAT_FILE"

    echo "====================================" >> "$LOG_FILE"
    echo >> "$LOG_FILE"
}

# Process ASV pipelines
for pipeline in ${PIPELINES["ASV"]}; do
    process_pipeline "$pipeline" "ASV" ""
done

# Process OTU pipelines
for pipeline in ${PIPELINES["OTU"]}; do
    process_pipeline "$pipeline" "OTU" "_otu"
done
