#!/bin/bash
#SBATCH -p rome
#SBATCH -t 12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/log/taxonomy/IDTAXA_output_%j.out

# Classifies a pipeline’s ASVs with IDTAXA against the V4-trained SILVA and Greengenes2 databases.

############################################################################################
#                              Pipeline configuration                                      #
############################################################################################

IDTAXA="/gpfs/work2/0/gs/danique/pipelines/scripts/classify/IDTAXA.R"
SILVA="https://www2.decipher.codes/data/Downloads/TrainingSets/SILVA_SSU_r138_2_2024.RData"
SILVA_NAME=$(basename $SILVA)

V4="/gpfs/work2/0/gs/danique/pipelines/simulations/input/db/trainset.RData"
GG_V4="/gpfs/work2/0/gs/danique/pipelines/ref/trainset_gg.RData"


############################################################################################
#                                   Run configuration                                      #
############################################################################################

# Assess number of assigned threads
if [ -z "$SLURM_CPUS_PER_TASK" ]; then
  THREADS=96
else
  # THREADS=$SLURM_CPUS_PER_TASK
  THREADS=$(($SLURM_NTASKS * $SLURM_CPUS_PER_TASK))
fi

# Access r libraries
R_LIBS=/projects/0/gs/danique/rpackages:$R_LIBS
export R_LIBS

# Load R module (pre-installed in Snellius)
module load 2024
module load R/4.4.2-gfbf-2024a

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

mkdir -p ${OUTPUT_DIR}
mkdir -p ${REF_DIR}

SEQ=$(find ${OUTPUT_DIR} -type f -name "$SEQNAME" | head -n 1)

if [ -z "$SEQ" ]; then
  echo "Error: Sequence file not found in ${OUTPUT_DIR}"
  exit 1
fi

echo "Dataset: $DATASET"
echo "Pipeline: $PIPELINE"
echo "Sequence file: $SEQ"
echo

# Install reference database if not already installed
# Silva (v138)
if [ ! -f "${REF_DIR}/${SILVA_NAME}" ]; then
  wget -P ${REF_DIR} ${SILVA}
fi

# # Run ID taxa R script
# echo "Run ID taxa with pretained silva trainingset"
# Rscript --verbose $IDTAXA -s $SEQ -t ${REF_DIR}/${SILVA_NAME} -o ${OUTPUT_DIR} -c ${THREADS} -n "pretrained"

# Run ID taxa R script
echo "Run ID taxa with V4 training set"
Rscript --verbose $IDTAXA -s $SEQ -t $V4 -o ${OUTPUT_DIR} -c ${THREADS} -n "V4"

# Run ID taxa R script wit greengenes
echo "Run ID taxa with GreenGenes V4 training set"
Rscript --verbose $IDTAXA -s $SEQ -t $GG_V4 -o ${OUTPUT_DIR} -c ${THREADS} -n "gg_V4"
