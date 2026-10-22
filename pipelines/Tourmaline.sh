#!/bin/bash
#SBATCH -p rome
#SBATCH -t 2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/pipelines/Tourmaline/log/Tourmaline_output_%j.out

############################################################################################
#                              Pipeline configuration                                      #
############################################################################################

PIPELINE="Tourmaline"
SILVASEQ="https://data.qiime2.org/2023.2/common/silva-138-99-seqs-515-806.qza"
SILVASEQ_NAME=$(basename ${SILVASEQ})
SILVATAX="https://data.qiime2.org/2023.2/common/silva-138-99-tax-515-806.qza"
SILVATAX_NAME=$(basename ${SILVATAX})

############################################################################################
#                                   Run configuration                                      #
############################################################################################

if [ -z "$SLURM_CPUS_PER_TASK" ]; then
  THREADS=96 
else
  # THREADS=$SLURM_CPUS_PER_TASK
  THREADS=$(($SLURM_NTASKS * $SLURM_CPUS_PER_TASK))
fi

CONDA_BASE=$(conda info --base)
CONDA_SH="$CONDA_BASE/etc/profile.d/conda.sh"
source $CONDA_SH

# Set variables (which data to use, where to find it, and which primers are attached)
# Set default values
DEFAULT_DATASET="mock5"
DEFAULT_DATA_DIR="/gpfs/work2/0/gs/danique/mockrobiota/mock5/fastq_demux/"
DEFAULT_FORWARD_PRIMER="GTGCCAGCMGCCGCGGTAA"
DEFAULT_REVERSE_PRIMER="GGACTACHVGGGTWTCTAAT"
DEFAULT_PRIMER_ATTACHED=0
DEFAULT_TRUNC_F=190
DEFAULT_TRUNC_R=135
DEFAULT_PATTERN="_L001__R[12]_001.fastq.gz"

while getopts ":d:i:f:r:a:t:u:e:h:" opt; do
  case $opt in
    d) DATASET="$OPTARG";;
    i) DATA_DIR="$OPTARG";;
    f) FORWARD_PRIMER="$OPTARG";;
    r) REVERSE_PRIMER="$OPTARG";;
    a) PRIMER_ATTACHED="$OPTARG";;
    t) TRUNC_F="$OPTARG";;
    u) TRUNC_R="$OPTARG";; 
    e) PATTERN="$OPTARG";;    
    h)
      echo "Usage: \$0 [options]"
      echo "Options:"
      echo "  -d <dataset>            Specify dataset. This will be the directory name in which intermediate (intermediate/<DATASET>) and output (output/<DATASET>) files are stored"
      echo "  -i <data_dir>           Specify data directory. This directory should contain all fastq.gz files to be processed"
      echo "  -f <forward_primer>     Specify forward primer"
      echo "  -r <reverse_primer>     Specify reverse primer"
      echo "  -a <primer_attached>    Specify whether primers are still attached."
      echo "  -t <truncate_forward>   Specify at which position to truncate forward reads" 
      echo "  -u <truncate_reverse>   Specify at which position to truncate reverse reads"  
      echo "  -e <pattern>            Specify the pattern of the fastq files"     
      echo "  -h                      Display this help message"
      exit 0
      ;;
    \?) echo "Invalid option: -$OPTARG"; exit 1;;
  esac
done

# Set variables to default values if not provided
DATASET=${DATASET:-$DEFAULT_DATASET}
DATA_DIR=${DATA_DIR:-$DEFAULT_DATA_DIR}
FORWARD_PRIMER=${FORWARD_PRIMER:-$DEFAULT_FORWARD_PRIMER}
REVERSE_PRIMER=${REVERSE_PRIMER:-$DEFAULT_REVERSE_PRIMER}
PRIMER_ATTACHED=${PRIMER_ATTACHED:-$DEFAULT_PRIMER_ATTACHED}
TRUNC_F=${TRUNC_F:-$DEFAULT_TRUNC_F}
TRUNC_R=${TRUNC_R:-$DEFAULT_TRUNC_R}
PATTERN=${PATTERN:-$DEFAULT_PATTERN}

# Pipeline Directories
BASE_DIR="/gpfs/work2/0/gs/danique/pipelines/pipelines/${PIPELINE}"
INTERMEDIATE_DIR="${BASE_DIR}/intermediate/${DATASET}"
OUTPUT_DIR="${BASE_DIR}/output/${DATASET}"
REF_DIR="${BASE_DIR}/intermediate"

mkdir -p ${INTERMEDIATE_DIR}
mkdir -p ${OUTPUT_DIR}

############################################################################################
#                                   Install databases                                      #
############################################################################################

if [ ! -f "${REF_DIR}/${SILVASEQ_NAME}" ]; then
  wget -P "${REF_DIR}/" ${SILVASEQ}
fi  
  
if [ ! -f "${REF_DIR}/${SILVATAX_NAME}" ]; then
  wget -P "${REF_DIR}/" ${SILVATAX}
fi    

############################################################################################
#                                            Link Data                                     #
############################################################################################

# Everytime you run Tourmaline on a new dataset/same dataset with new parameters, clone tourmaline from github, Tihs allows snakemake to run smoothly without having to account for previous inputs and output
cd $BASE_DIR 

if [ -d "${BASE_DIR}/tourmaline" ]; then
  rm -Rf "${BASE_DIR}/tourmaline" 
fi

if [ -d "${BASE_DIR}/${DATASET}" ]; then
  rm -Rf "${BASE_DIR}/${DATASET}" 
fi

git clone -b V1 https://github.com/aomlomics/tourmaline.git # V1 to ensure version stability 

# Rename tourmaline dir to the name of the project
mv tourmaline $DATASET
cd "${BASE_DIR}/${DATASET}"
      
# Link reference database
if [ ! -d  ./01-imported ]; then
  mkdir ./01-imported
fi

ln -sf "${REF_DIR}/${SILVASEQ_NAME}" ./01-imported/refseqs.qza
ln -sf "${REF_DIR}/${SILVATAX_NAME}" ./01-imported/reftax.qza

# empty the data directory (remove test data)
if [ -d  ./00-data/ ]; then
  if [ "$(ls -A ./00-data/)" ]; then
    rm -rf ./00-data/* 
  fi  
fi

mkdir ./00-data/fastq

# Link fastq files
find "$DATA_DIR" -maxdepth 1 -name "*.fastq.gz" -exec ln -s {} ./00-data/fastq \;

############################################################################################
#                                Create manifest_file                                      #
############################################################################################

# To ensure sample IDs are properly parsed 
sed -i "s/('_R1_' in x) or ('_R2_' in x)/('_R1' in x) or ('_R2' in x)/g" ./scripts/create_manifest_from_fastq_directory.py  
sed -i "s/_S\[0-9\]{1,3}_L001_R\[12\]_001.fastq.gz/${PATTERN}/g" ./scripts/create_manifest_from_fastq_directory.py

# Create manifest file
./scripts/create_manifest_from_fastq_directory.py "${BASE_DIR}/${DATASET}/00-data/fastq" ./00-data/manifest_pe.csv ./00-data/manifest_se.csv

# Create metadata file 
awk -F, -v OFS="\t" '{ $1=$1; print }' 00-data/manifest_se.csv | \
awk -F'\t' '{print $1}' | \
awk -F, -v OFS="\t" 'BEGIN {print "sample-id", "composition"} {if (NR > 1) {if ($1 ~ /Even/) {print $0, "even"} else if ($1 ~ /Staggered/) {print $0, "staggered"} else {print $0, "unknown"}}}' > 00-data/metadata.tsv

############################################################################################
#                                 Adjust configuration file                                #
############################################################################################

# Adjust configuration file for dataset
CONFIG_FILE="${BASE_DIR}/${DATASET}/config.yaml"

# Data specification
yq --arg TRUNC_F "$TRUNC_F" '.dada2pe_trunc_len_f = ($TRUNC_F | tonumber)' -i -y ${CONFIG_FILE} # trim forward
yq --arg TRUNC_R "$TRUNC_R" '.dada2pe_trunc_len_r = ($TRUNC_R | tonumber)' -i -y ${CONFIG_FILE} # trim reverse
yq --arg THREADS "$THREADS" '.dada2pe_threads = ($THREADS | tonumber)' -i -y ${CONFIG_FILE}
yq --arg THREADS "$THREADS" '.dada2se_threads = ($THREADS | tonumber)' -i -y ${CONFIG_FILE}
yq --arg THREADS "$THREADS" '.deblur_threads = ($THREADS | tonumber)' -i -y ${CONFIG_FILE}
yq --arg THREADS "$THREADS" '.alignment_threads = ($THREADS | tonumber)' -i -y ${CONFIG_FILE}
yq --arg THREADS "$THREADS" '.feature_classifier_threads = ($THREADS | tonumber)' -i -y ${CONFIG_FILE}
yq --arg THREADS "$THREADS" '.phylogeny_fasttree_threads = ($THREADS | tonumber)' -i -y ${CONFIG_FILE}
yq --arg THREADS "$THREADS" '.diversity_core_metrics_phylogenetic_threads = ($THREADS | tonumber)' -i -y ${CONFIG_FILE}
yq --arg THREADS "$THREADS" '.other_threads = ($THREADS | tonumber)' -i -y ${CONFIG_FILE}

yq '.beta_group_column = "composition"' -i -y ${CONFIG_FILE}

if [ "$PRIMER_ATTACHED" -eq 1 ]; then
  yq --arg LEN_F "${#FORWARD_PRIMER}" '.dada2pe_trim_left_f = ($LEN_F | tonumber)' -i -y ${CONFIG_FILE} 
  yq --arg LEN_R "${#REVERSE_PRIMER}" '.dada2pe_trim_left_r = ($LEN_R | tonumber)' -i -y ${CONFIG_FILE} 
fi

############################################################################################
#                                      Run Pipeline                                        #
############################################################################################

conda activate snakemake

snakemake --use-conda dada2_pe_denoise
snakemake --use-conda dada2_pe_taxonomy_unfiltered
snakemake --use-conda dada2_pe_diversity_unfiltered
snakemake --use-conda dada2_pe_report_unfiltered

# Filtering (default filters)

if [ -f "./02-output-dada2-pe-unfiltered/02-alignment-tree/repseqs_to_filter_outliers.tsv" ]; then
  cp ./02-output-dada2-pe-unfiltered/02-alignment-tree/repseqs_to_filter_outliers.tsv ./00-data/repseqs_to_filter_dada2-pe.tsv
else 
  echo -e "featureid\toutlier" > ./00-data/repseqs_to_filter_dada2-pe.tsv
  echo -e "abcd1234\ttrue" >> ./00-data/repseqs_to_filter_dada2-pe.tsv # add nonestisting feature id
fi

echo "#SampleID" > ./00-data/samples_to_filter_dada2-pe.tsv
echo "TMP123" >> ./00-data/samples_to_filter_dada2-pe.tsv # Add a nonexisting sample id

snakemake --use-conda dada2_pe_taxonomy_filtered 
snakemake --use-conda dada2_pe_diversity_filtered
snakemake --use-conda dada2_pe_report_filtered 

############################################################################################
#                                      Export output                                       #
############################################################################################

# Unfiltered
awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' ./02-output-dada2-pe-unfiltered/00-table-repseqs/repseqs.fasta > "${OUTPUT_DIR}/dna-sequences.fasta"
cp ./02-output-dada2-pe-unfiltered/01-taxonomy/taxonomy.tsv "${OUTPUT_DIR}/taxonomy.tsv"
cp ./02-output-dada2-pe-unfiltered/01-taxonomy/asv_taxa_sample_table.tsv "${OUTPUT_DIR}/feature-table.tsv"

# Filtered
awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' ./02-output-dada2-pe-filtered/00-table-repseqs/repseqs.fasta > "${OUTPUT_DIR}/dna-sequences_filt.fasta"
cp ./02-output-dada2-pe-unfiltered/01-taxonomy/taxonomy.tsv "${OUTPUT_DIR}/taxonomy_filt.tsv"
cp ./02-output-dada2-pe-unfiltered/01-taxonomy/asv_taxa_sample_table.tsv "${OUTPUT_DIR}/feature-table_filt.tsv"


# Link tourmaline dir to intermediate dir
if [ -L "${INTERMEDIATE_DIR}/${DATASET}" ]; then
  rm "${INTERMEDIATE_DIR}/${DATASET}"
fi

ln -s "${BASE_DIR}/${DATASET}" "${INTERMEDIATE_DIR}/${DATASET}"

############################################################################################
#                                      Reformat output                                     #
############################################################################################

# dna-sequences
cp "${OUTPUT_DIR}/dna-sequences.fasta" "${OUTPUT_DIR}/$(basename dna-sequences.fasta .fasta)_post.fasta"

# feature table 
cut -f1,5- "${OUTPUT_DIR}/feature-table.tsv" | \
sed '1s/^featureid/ID/' > "${OUTPUT_DIR}/$(basename feature-table.tsv .tsv)_post.tsv"

# taxonomy
awk '{ NF--; print }' OFS='\t' "${OUTPUT_DIR}/taxonomy.tsv" | \
sed 's/;/\t/g;s/.\?__//g' | sed '1d' | \
awk 'BEGIN {OFS="\t"; print "ID", "kingdom", "phylum", "class", "order", "family", "genus", "species"} 
     {print $1, "NA", $2, $3, $4, $5, $6, $7}' > "${OUTPUT_DIR}/$(basename taxonomy.tsv .tsv)_post.tsv" 


############################################################################################
#                                      Descriptive stats                                   #
############################################################################################

# Write basic descriptive stats to txt file
numasv=$(awk '/^>/' "${OUTPUT_DIR}/dna-sequences_post.fasta" | wc -l) 
numgenus=$(awk -F'\t' '{print $7 }' "${OUTPUT_DIR}/taxonomy_post.tsv" | sort | uniq -c | wc -l)

echo "asv: ${numasv}" > "${OUTPUT_DIR}/stat.txt"
echo "genera: ${numgenus}" >> "${OUTPUT_DIR}/stat.txt"

awk -F'\t' '{print $7 }' "${OUTPUT_DIR}/taxonomy_post.tsv" | sort | uniq -c | while read line; do
  echo "$line" >> "${OUTPUT_DIR}/stat.txt"
done