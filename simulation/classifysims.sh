#!/bin/bash
#SBATCH -p rome
#SBATCH -t 05:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/log/taxonomy/classify_mock_%j.out

# Classifies the representative sequences (to simulate) with Q2-NBC and IDTAXA.

while getopts ":c:d:" opt; do
  case $opt in
    c) CLASSIFIER="$OPTARG";;
    d) DATASET="$OPTARG";;
    \?) echo "Invalid option: -$OPTARG"; exit 1;;
  esac
done

echo "Running ${CLASSIFIER}..." 
echo "Dataset ${DATASET}..." 

if [ -z "$SLURM_CPUS_PER_TASK" ]; then
  THREADS=24
else
  # THREADS=$SLURM_CPUS_PER_TASK
  THREADS=$(($SLURM_CPUS_PER_TASK))
fi

OUTPUT_DIR="/projects/0/gs/danique/pipelines/mock/tax"
REF_DIR="/gpfs/work2/0/gs/danique/pipelines/ref"
mkdir -p $OUTPUT_DIR

###################################
## QIIME2 Naive Bayes Classifier ##
###################################

if grep -q "q2" <<< $CLASSIFIER ; then

    SILVA="https://data.qiime2.org/2024.2/common/silva-138-99-515-806-nb-classifier.qza"
    SILVA_NAME=$(basename $SILVA)

    if [ ! -f "${REF_DIR}/${SILVA_NAME}" ]; then
        wget -P ${REF_DIR} ${SILVA}
    fi

    # # load QIIME2 software (pre-installed in snellius)
    # module load 2023
    # module load QIIME2/2023.7

    CONDA_BASE=$(conda info --base)
    CONDA_SH="$CONDA_BASE/etc/profile.d/conda.sh"
    source $CONDA_SH

    for m in $DATASET; do

        if grep -q "mock" <<< $m ; then
            SEQ="/projects/0/gs/danique/mockrobiota/${m}/info/expected-sequences_num.fasta"
            sed 's/ /_/g' $SEQ > ${OUTPUT_DIR}/tmp.fasta
        else
            SEQ="/projects/0/gs/danique/pipelines/simulations/input/tosimulate_resolvable.fasta"
            sed 's/ /_/g' $SEQ > ${OUTPUT_DIR}/tmp.fasta
        fi

        SEQ="${OUTPUT_DIR}/tmp.fasta"

        conda activate qiime2-2024.2

        # Import taxonomy tsv into qiime2 artifact
        qiime tools import \
        --input-path $SEQ \
        --output-path "${OUTPUT_DIR}/rep-seqs_${m}.qza" \
        --type 'FeatureData[Sequence]'

        echo "Classify sequences for ${m} using QIIME2..."
        # Silva
        qiime feature-classifier classify-sklearn \
        --i-classifier "${REF_DIR}/${SILVA_NAME}" \
        --i-reads "${OUTPUT_DIR}/rep-seqs_${m}.qza" \
        --o-classification "${OUTPUT_DIR}/q2_${m}.qza" \
        --p-n-jobs ${THREADS}

        echo "Export..."

        qiime tools export \
        --input-path "${OUTPUT_DIR}/q2_${m}.qza" \
        --output-path  "${OUTPUT_DIR}/"

        # Reformat
        awk -F'\t' '{ $NF=""; print $0 }' OFS='\t' ${OUTPUT_DIR}/taxonomy.tsv | \
        sed 's/;/\t/g;s/.\?__//g' | \
        awk 'BEGIN {OFS="\t"; print "ID", "kingdom", "phylum", "class", "order", "family", "genus", "species"} {if (NR > 1) print}' > "${OUTPUT_DIR}/q2_${m}.tsv"

        rm "${OUTPUT_DIR}/taxonomy.tsv"
    done
fi

###################################
##            IDTaxa             ##
###################################

if  grep -q "idtaxa" <<< $CLASSIFIER ; then 

    IDTAXA="/gpfs/work2/0/gs/danique/pipelines/scripts/classify/IDTAXA.R"
    V4="/gpfs/work2/0/gs/danique/pipelines/simulations/input/db/trainset.RData"

    # Access r libraries
    R_LIBS=/projects/0/gs/danique/rpackages:$R_LIBS
    export R_LIBS

    # Load R module (pre-installed in Snellius)
    module load 2024
    module load R/4.4.2-gfbf-2024a

    for m in $DATASET; do

        if grep -q "mock" <<< $m ; then
            SEQ="/projects/0/gs/danique/mockrobiota/${m}/info/expected-sequences_num.fasta"
            sed 's/ /_/g' $SEQ > ${OUTPUT_DIR}/tmp.fasta
        else
            SEQ="/projects/0/gs/danique/pipelines/simulations/input/tosimulate_resolvable.fasta"
            sed 's/ /_/g' $SEQ > ${OUTPUT_DIR}/tmp.fasta
        fi

        # Run ID taxa R script
        Rscript --verbose $IDTAXA -s $SEQ -t ${V4} -o ${OUTPUT_DIR} -c ${THREADS} -n "v4"
        mv ${OUTPUT_DIR}/idtaxa_v4.tsv ${OUTPUT_DIR}/idtaxa_${m}_v4.tsv    
        mv ${OUTPUT_DIR}/idtaxa_confidence_v4.tsv ${OUTPUT_DIR}/idtaxa_confidence_${m}_v4.tsv  

    done
fi

