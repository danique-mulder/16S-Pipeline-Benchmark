#!/bin/bash
#SBATCH -p rome
#SBATCH -t 2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/pipelines/APSCALE/log/APSCALE_output_%j.out

############################################################################################
#                              Pipeline configuration                                      #
############################################################################################

PIPELINE="APSCALE"

############################################################################################
#                                   Run configuration                                      #
############################################################################################

if [ -z "$SLURM_CPUS_PER_TASK" ]; then
  THREADS=96
else
  # THREADS=$SLURM_CPUS_PER_TASK
  THREADS=$(($SLURM_NTASKS * $SLURM_CPUS_PER_TASK))
fi
# test

# source vsearch
export PATH=$PATH:/gpfs/work2/0/gs/danique/pipelines/pipelines/VSEARCH/vsearch-2.28.1/bin/

# Set variables (which data to use, where to find it, and which primers are attached)
# Set default values
DEFAULT_DATASET="mock5"
DEFAULT_DATA_DIR="/gpfs/work2/0/gs/danique/mockrobiota/mock5/fastq_demux/primer"
DEFAULT_FORWARD_PRIMER="GTGCCAGCMGCCGCGGTAA"
DEFAULT_REVERSE_PRIMER="GGACTACHVGGGTWTCTAAT"
DEFAULT_PRIMER_ATTACHED=1
DEFAULT_REMOVEOLD=0

while getopts ":d:i:f:r:a:y:h:" opt; do
  case $opt in
    d) DATASET="$OPTARG";;
    i) DATA_DIR="$OPTARG";;
    f) FORWARD_PRIMER="$OPTARG";;
    r) REVERSE_PRIMER="$OPTARG";;
    a) PRIMER_ATTACHED="$OPTARG";;
    y) REMOVEOLD="$OPTARG";;
    h)
      echo "Usage: \$0 [options]"
      echo "Options:"
      echo "  -d <dataset>            Specify dataset. This will be the directory name in which intermediate (intermediate/<DATASET>) and output (output/<DATASET>) files are stored"
      echo "  -i <data_dir>           Specify data directory. This directory should contain all fastq.gz files to be processed"
      echo "  -f <forward_primer>     Specify forward primer"
      echo "  -r <reverse_primer>     Specify reverse primer"
      echo "  -a <primer_attached>    Specify whether primers are still attached."
      echo "  -y <remove_old>         Remove old files (1) or keep them (0). Default is 0"
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
REMOVEOLD=${REMOVEOLD:-$DEFAULT_REMOVEOLD}

if [ "$PRIMER_ATTACHED" -eq 0 ]; then
  echo "ERROR - Sequences need to have primers attached"
  exit 1
fi 

# Pipeline Directories
BASE_DIR="/gpfs/work2/0/gs/danique/pipelines/pipelines/${PIPELINE}"
INTERMEDIATE_DIR="${BASE_DIR}/intermediate/${DATASET}"
OUTPUT_DIR="${BASE_DIR}/output/${DATASET}"

mkdir -p ${INTERMEDIATE_DIR}
mkdir -p ${OUTPUT_DIR}

# Remove files from previous run if REMOVEOLD is set to 1
if [ $REMOVEOLD -eq 1 ]; then
  echo "removing intermediate files from previous run on ${DATASET}"
  rm -rf "${INTERMEDIATE_DIR}"/* # ensure dir is empty
fi

############################################################################################
#                                Adjust configuration                                      #
############################################################################################

#NOTE: the xlsx configuration file needs to be checked manually.

# Create project
cd $INTERMEDIATE_DIR

if [ -d "$INTERMEDIATE_DIR/${DATASET}_apscale" ]; then
  rm -rf "$INTERMEDIATE_DIR/${DATASET}_apscale"
fi

apscale --create_project ${DATASET}
PROJECT_DIR="${INTERMEDIATE_DIR}/${DATASET}_apscale"
CONFIG="${PROJECT_DIR}/Settings_${DATASET}.xlsx"

# Add symbolic link to data
find "$DATA_DIR" -maxdepth 1 -name "*.fastq.gz" -exec ln -s {} "${PROJECT_DIR}/2_demultiplexing/data" \;

# Copy the configuration file (check/adjust manually)
cp "${BASE_DIR}/Settings_mock_primer.xlsx" "${PROJECT_DIR}/Settings_${DATASET}.xlsx"

############################################################################################
#                                      Run Pipeline                                        #
############################################################################################

# Run pipeline
apscale --run_apscale $PROJECT_DIR

############################################################################################
#                                      Export output                                       #
############################################################################################

# dna-sequences
awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' "${PROJECT_DIR}/8_esv_table/${DATASET}_ESVs.fasta" > "${OUTPUT_DIR}/dna-sequences.fasta"

# feature table
in2csv "${PROJECT_DIR}/8_esv_table/${DATASET}_ESV_table.xlsx" | csvformat -T > "${OUTPUT_DIR}/feature-table.tsv"

############################################################################################
#                                      Reformat output                                     #
############################################################################################

# dna-sequences
cp "${OUTPUT_DIR}/dna-sequences.fasta" "${OUTPUT_DIR}/$(basename dna-sequences.fasta .fasta)_post.fasta"

# feature table
awk -F'\t' -v OFS='\t' '{for (i=2; i<NF; i++) printf "%s%s", $i, (i<NF-1?OFS:"\n")}' "${OUTPUT_DIR}/feature-table.tsv" | \
sed '1s/^unique_//' | \
awk 'BEGIN {FS=OFS="\t"} 
     NR==1 {for (i=1; i<=NF; i++) {split($i, a, "_"); $i=a[1]}} 
     {print}' > "${OUTPUT_DIR}/$(basename feature-table.tsv .tsv)_post.tsv"


############################################################################################
#                                      Descriptive stats                                   #
############################################################################################

## ASV
# Write basic descriptive stats to txt file
numasv=$(awk '/^>/' "${OUTPUT_DIR}/dna-sequences_post.fasta" | wc -l) 
echo "asv: ${numasv}" > "${OUTPUT_DIR}/stat.txt"