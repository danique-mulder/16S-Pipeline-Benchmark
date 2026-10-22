#!/bin/bash
#SBATCH -p rome
#SBATCH -t 1-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/log/taxonomy/q2_naive_output_%j.out

############################################################################################
#                              Pipeline configuration                                      #
############################################################################################

SILVA="https://data.qiime2.org/2024.2/common/silva-138-99-515-806-nb-classifier.qza"
SILVA_NAME=$(basename $SILVA)
GG="https://data.qiime2.org/classifiers/greengenes/gg_2022_10_backbone.v4.nb.qza"
GG_NAME=$(basename $GG)

############################################################################################
#                                   Run configuration                                      #
############################################################################################

# load QIIME2 software (pre-installed in snellius)
# module load 2023
# module load QIIME2/2023.7

CONDA_BASE=$(conda info --base)
CONDA_SH="$CONDA_BASE/etc/profile.d/conda.sh"
source $CONDA_SH


# Set values for variables based on input
while getopts ":d:p:s:" opt; do
  case $opt in
    d) DATASET="$OPTARG";;
    p) PIPELINE="$OPTARG";;
    s) SEQNAME="$OPTARG";;
    h)
      echo "Usage: \$0 [options]"
      echo "Options:"
      echo "  -d <dataset>          Specify from which dataset the seqeunces are derived"
      echo "  -p <pipeline>         Specify from which pipeline the sequences are derived"
      echo "  -s <sequence_file>    Specify the sequence file path (fasta) to use"
      echo "  -h                    Display this help message"
      exit 0
      ;;
    \?) echo "Invalid option: -$OPTARG"; exit 1;;
  esac
done

# Install reference database if not already installed
OUTPUT_DIR="/gpfs/work2/0/gs/danique/pipelines/pipelines/${PIPELINE}/output/${DATASET}"
REF_DIR="/gpfs/work2/0/gs/danique/pipelines/ref"
SEQ=$(find ${OUTPUT_DIR} -type f -name "$SEQNAME" | head -n 1)
mkdir -p ${REF_DIR}

if [ -z "$SEQ" ]; then
  echo "Error: Sequence file not found in ${OUTPUT_DIR}"
  exit 1
fi

# Sometimes sequence files contain lowercase letters, which can cause issues with QIIME2 import
tmpfile=$(mktemp /tmp/tempfile_XXXXXX.fasta)
awk '{if($0 !~ /^>/) print toupper($0); else print $0}' $SEQ > ${tmpfile}


if [[ "$SEQ" == *_otu* ]]; then
  SUFFIX="_otu"
  TYPE="OTU"
else
  SUFFIX=""
  TYPE="ASV"
fi

if [[ "$SEQ" == *concat* ]]; then
  CONCAT="_concat"
else
  CONCAT=""
fi

mkdir -p ${OUTPUT_DIR}/tmp${CONCAT}${SUFFIX}

echo "Dataset: $DATASET"
echo "Pipeline: $PIPELINE"
echo "Sequence file: $SEQ"
echo

# Install reference database if not already installed
# Silva (v138)
if [ ! -f "${REF_DIR}/${SILVA_NAME}" ]; then
  wget -P ${REF_DIR} ${SILVA}
fi

# GreenGenes2
if [ ! -f "${REF_DIR}/${GG_NAME}" ]; then
  wget -P ${REF_DIR} ${GG}
fi

# Load env
conda activate qiime2-2024.2

# Import taxonomy tsv into qiime2 artifact
qiime tools import \
  --input-path ${tmpfile} \
  --output-path "${OUTPUT_DIR}/tmp${CONCAT}${SUFFIX}/rep-seqs${CONCAT}${SUFFIX}.qza" \
  --type 'FeatureData[Sequence]'

rm ${tmpfile}

qiime feature-table tabulate-seqs \
  --i-data "${OUTPUT_DIR}/tmp${CONCAT}${SUFFIX}/rep-seqs${CONCAT}${SUFFIX}.qza" \
  --o-visualization "${OUTPUT_DIR}/rep-seqs${CONCAT}${SUFFIX}"

############################################################################################
#                                  Run Silva Classifier                                    #
############################################################################################

# # Silva
# qiime feature-classifier classify-sklearn \
#   --i-classifier "${REF_DIR}/${SILVA_NAME}" \
#   --i-reads "${OUTPUT_DIR}/tmp${CONCAT}${SUFFIX}/rep-seqs${CONCAT}${SUFFIX}.qza" \
#   --o-classification "${OUTPUT_DIR}/tmp${CONCAT}${SUFFIX}/q2_nb_silva${CONCAT}${SUFFIX}.qza" \
#   --p-n-jobs 1 

# # To avoid two processes writing the same file at the same time 
# MY_TMPDIR="${TMPDIR}/q2_nb_silva${CONCAT}${SUFFIX}/${PIPELINE}/"
# mkdir -p "$MY_TMPDIR"

# qiime tools export \
#   --input-path "${OUTPUT_DIR}/tmp${CONCAT}${SUFFIX}/q2_nb_silva${CONCAT}${SUFFIX}.qza" \
#   --output-path  "${MY_TMPDIR}/"

# mv "${MY_TMPDIR}/taxonomy.tsv" "${OUTPUT_DIR}/tmp${CONCAT}${SUFFIX}/q2_nb_silva${CONCAT}${SUFFIX}.tsv"
# rm -R "$MY_TMPDIR"


# ############################################################################################
# #                                 Export taxonomic asignments                              #
# ############################################################################################

# if [ -f "${OUTPUT_DIR}/tmp${CONCAT}${SUFFIX}/q2_nb_silva${CONCAT}${SUFFIX}.tsv" ]; then
#   # Silva
#   awk -F'\t' '{ $NF=""; print $0 }' OFS='\t' "${OUTPUT_DIR}/tmp${CONCAT}${SUFFIX}/q2_nb_silva${CONCAT}${SUFFIX}.tsv" | \
#   sed 's/;/\t/g;s/.\?__//g' | \
#   awk 'BEGIN {OFS="\t"; print "ID", "kingdom", "phylum", "class", "order", "family", "genus", "species"} {if (NR > 1) print}' > "${OUTPUT_DIR}/q2_nb_silva${CONCAT}${SUFFIX}.tsv"

#   echo
#   head "${OUTPUT_DIR}/q2_nb_silva${CONCAT}${SUFFIX}.tsv"
#   echo

#   # Extract the first and last columns and replace headers
#   awk -F'\t' 'NR==1 {OFS="\t"; print "ID", "confidence"; next} {print $1, $NF}' "${OUTPUT_DIR}/tmp${CONCAT}${SUFFIX}/q2_nb_silva${CONCAT}${SUFFIX}.tsv" > "${OUTPUT_DIR}/q2_nb_silva_confidence${CONCAT}${SUFFIX}.tsv"

#   echo
#   head "${OUTPUT_DIR}/q2_nb_silva_confidence${CONCAT}${SUFFIX}.tsv"
#   echo

# else
#   echo "Error: Silva classification file not found"
#   exit 1
# fi

############################################################################################
#                               Run GreenGenes Classifier                                  #
############################################################################################

# gg
qiime feature-classifier classify-sklearn \
  --i-classifier "${REF_DIR}/${GG_NAME}" \
  --i-reads "${OUTPUT_DIR}/tmp${CONCAT}${SUFFIX}/rep-seqs${CONCAT}${SUFFIX}.qza" \
  --o-classification "${OUTPUT_DIR}/tmp${CONCAT}${SUFFIX}/q2_nb_gg${CONCAT}${SUFFIX}.qza" \
  --p-n-jobs 1 

# To avoid two processes writing the same file at the same time 
MY_TMPDIR="${TMPDIR}/q2_nb_gg${CONCAT}${SUFFIX}/${PIPELINE}/"
mkdir -p "$MY_TMPDIR"

qiime tools export \
  --input-path "${OUTPUT_DIR}/tmp${CONCAT}${SUFFIX}/q2_nb_gg${CONCAT}${SUFFIX}.qza" \
  --output-path  "${MY_TMPDIR}/"

mv "${MY_TMPDIR}/taxonomy.tsv" "${OUTPUT_DIR}/tmp${CONCAT}${SUFFIX}/q2_nb_gg${CONCAT}${SUFFIX}.tsv"
rm -R "$MY_TMPDIR"


############################################################################################
#                                 Export taxonomic asignments                              #
############################################################################################

if [ -f "${OUTPUT_DIR}/tmp${CONCAT}${SUFFIX}/q2_nb_gg${CONCAT}${SUFFIX}.tsv" ]; then
  # gg
  awk -F'\t' '{ $NF=""; print $0 }' OFS='\t' "${OUTPUT_DIR}/tmp${CONCAT}${SUFFIX}/q2_nb_gg${CONCAT}${SUFFIX}.tsv" | \
  sed 's/;/\t/g;s/.\?__//g' | \
  awk 'BEGIN {OFS="\t"; print "ID", "kingdom", "phylum", "class", "order", "family", "genus", "species"} {if (NR > 1) print}' > "${OUTPUT_DIR}/q2_nb_gg${CONCAT}${SUFFIX}.tsv"

  echo
  head "${OUTPUT_DIR}/q2_nb_gg${CONCAT}${SUFFIX}.tsv"
  echo

  # Extract the first and last columns and replace headers
  awk -F'\t' 'NR==1 {OFS="\t"; print "ID", "confidence"; next} {print $1, $NF}' "${OUTPUT_DIR}/tmp${CONCAT}${SUFFIX}/q2_nb_gg${CONCAT}${SUFFIX}.tsv" > "${OUTPUT_DIR}/q2_nb_gg_confidence${CONCAT}${SUFFIX}.tsv"

  echo
  head "${OUTPUT_DIR}/q2_nb_gg_confidence${CONCAT}${SUFFIX}.tsv"
  echo

else
  echo "Error: GreenGenes2 classification file not found"
  exit 1
fi





# Remove intermediate files
rm -rf "${OUTPUT_DIR}/tmp${SUFFIX}"

conda deactivate
