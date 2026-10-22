#!/bin/bash
#SBATCH -p rome
#SBATCH -t 2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/pipelines/DADA2/log/DADA2_output_%j.out

############################################################################################
#                              Pipeline configuration                                      #
############################################################################################

PIPELINE="DADA2"
SILVA="https://zenodo.org/records/4587955/files/silva_nr99_v138.1_train_set.fa.gz"
SILVANAME=$(basename ${SILVA})
DADA=/gpfs/work2/0/gs/danique/pipelines/scripts/pipelines/DADA2.R

############################################################################################
#                                   Run configuration                                      #
############################################################################################

if [ -z "$SLURM_CPUS_PER_TASK" ]; then
  THREADS=96
else
  # THREADS=$SLURM_CPUS_PER_TASK
  THREADS=$(($SLURM_NTASKS * $SLURM_CPUS_PER_TASK))
fi

# # for cutadapt
# export PATH=$PATH://gpfs/home4/dmulder/miniforge3/bin/

# Access r libraries
R_LIBS=$HOME/rpackages:$R_LIBS
export R_LIBS

# Load R module (pre-installed in Snellius)
module load 2023
module load R/4.3.2-gfbf-2023a

# Set variables (which data to use, where to find it, and which primers are attached)
# Set default values
DEFAULT_DATASET="gemsim"
DEFAULT_DATA_DIR="/gpfs/work2/0/gs/danique/pipelines/data/gemsim/seq"
DEFAULT_FORWARD_PRIMER="GTGCCAGCMGCCGCGGTAA"
DEFAULT_REVERSE_PRIMER="GGACTACHVGGGTWTCTAAT"
DEFAULT_PRIMER_ATTACHED=0
DEFAULT_TRUNC_F=230
DEFAULT_TRUNC_R=230
DEFAULT_VERBOSE=1
DEFAULT_REMOVEOLD=0

while getopts ":d:i:f:r:a:t:u:v:y:h:" opt; do
  case $opt in
    d) DATASET="$OPTARG";;
    i) DATA_DIR="$OPTARG";;
    f) FORWARD_PRIMER="$OPTARG";;
    r) REVERSE_PRIMER="$OPTARG";;
    a) PRIMER_ATTACHED="$OPTARG";;
    t) TRUNC_F="$OPTARG";;
    u) TRUNC_R="$OPTARG";;    
    v) VERBOSE="$OPTARG";;
    y) REMOVEOLD="$OPTARG";;
    h)
      echo "Usage: \$0 [options]"
      echo "Options:"
      echo "  -d <dataset>          Specify dataset. This will be the directory name in which intermediate (intermediate/<DATASET>) and output (output/<DATASET>) files are stored"
      echo "  -i <data_dir>         Specify data directory. This directory should contain all fastq.gz files to be processed"
      echo "  -f <forward_primer>   Specify forward primer sequence"
      echo "  -r <reverse_primer>   Specify reverse primer sequence"
      echo "  -a <primer_attached>  Specify whether primers are still attached. NOTE: DADA2 requires sequences without adapters and primers."
      echo "  -t <truncate_forward> Specify at which position to truncate forward reads" 
      echo "  -u <truncate_reverse> Specify at which position to truncate reverse reads"       
      echo "  -v <verbose>          Print progress messages (0/1)"  
      echo "  -y <remove_old>       Remove old files (0/1)"
      echo "  -h                    Display this help message"
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
VERBOSE=${VERBOSE:-$DEFAULT_VERBOSE}
REMOVEOLD=${REMOVEOLD:-$DEFAULT_REMOVEOLD}

echo $TRUNC_F
echo $TRUNC_R

# Pipeline Directories
BASE_DIR="/gpfs/work2/0/gs/danique/pipelines/pipelines/${PIPELINE}"
INTERMEDIATE_DIR="${BASE_DIR}/intermediate/${DATASET}"
OUTPUT_DIR="${BASE_DIR}/output/${DATASET}"
REF_DIR="${BASE_DIR}/intermediate"

mkdir -p "${INTERMEDIATE_DIR}"
mkdir -p "${OUTPUT_DIR}"

# Remove files from previous run if REMOVEOLD is set to 1
if [ $REMOVEOLD -eq 1 ]; then
  echo "removing intermediate files from previous run on ${DATASET}"
  rm -rf "${INTERMEDIATE_DIR}/"*
fi

if [ $VERBOSE -eq 1 ]; then
  VERBOSE=TRUE
else
  VERBOSE=FALSE
fi

# Download silva DB, if not yet installed
if [ ! -f "${REF_DIR}/${SILVANAME}" ]; then
  wget -P $REF_DIR $SILVA
fi

ln -sf "${REF_DIR}/${SILVANAME}" ${INTERMEDIATE_DIR} # symbolic link to reference db

############################################################################################
#                                   Remove primers (if attached)                           #
############################################################################################

# Run cutadapt if primers are still attached (and if (part of) the files have not been trimmed)
if [ $PRIMER_ATTACHED -eq 1 ] && { [ ! -d "${INTERMEDIATE_DIR}/trimmed" ] || [ $(ls -1q ${INTERMEDIATE_DIR}/trimmed | wc -l) -ne $(ls -1q ${DATA_DIR}/*_R1*.fastq.gz | wc -l) ]; }; then
  # Create the trimmed directory if it doesn't exist
  mkdir -p "${INTERMEDIATE_DIR}/trimmed"
  mkdir -p "${INTERMEDIATE_DIR}/untrimmed"

  # Ensure the directories are empty  
  rm "${INTERMEDIATE_DIR}/trimmed"/* 
  rm "${INTERMEDIATE_DIR}/untrimmed"/* 

  # Loop through all R1 files
  for ffile in ${DATA_DIR}/*_R1*.fastq.gz; do
    rfile="${ffile/_R1/_R2}" # R2 file
    
    fout="${INTERMEDIATE_DIR}/trimmed/$(basename "${ffile}" .fastq.gz)_trimmed.fastq.gz"
    rout="${INTERMEDIATE_DIR}/trimmed/$(basename "${rfile}" .fastq.gz)_trimmed.fastq.gz"

    fout_untrimmed="${INTERMEDIATE_DIR}/untrimmed/$(basename "${ffile}" .fastq.gz)_untrimmed.fastq.gz"
    rout_untrimmed="${INTERMEDIATE_DIR}/untrimmed/$(basename "${rfile}" .fastq.gz)_untrimmed.fastq.gz"

    echo "R1: $(basename $ffile)        R2: $(basename $rfile)"
    cutadapt -g "^${FORWARD_PRIMER}" \
      -G "^${REVERSE_PRIMER}" \
      -o "${fout}" \
      -p "${rout}" \
      --untrimmed-output "${fout_untrimmed}" \
      --untrimmed-paired-output "${rout_untrimmed}" \
      -j "${THREADS}" \
      "${ffile}" "${rfile}" 
  done
  
  DATA_DIR="${INTERMEDIATE_DIR}/trimmed"
  
fi

# Find filename pattern for parsing of sample names
PATTERN=$(ls ${DATA_DIR} | grep -m 1 '_R1' | sed 's/.*R1/R1/')  # Get the first file with 'R1' in its name

############################################################################################
#                                      Run Pipeline                                        #
############################################################################################

# Run DADA2 R script
Rscript --verbose $DADA -p "${DATA_DIR}" -d "${INTERMEDIATE_DIR}" -o "${OUTPUT_DIR}" -c "${THREADS}" -a "${PATTERN}" -t "${TRUNC_F}" -u "${TRUNC_R}" -v "${VERBOSE}"

echo "path <- ${DATA_DIR}"
echo "tmp_path <- ${INTERMEDIATE_DIR}"
echo "output <- ${OUTPUT_DIR}"
echo "verbose <- ${VERBOSE}"
echo "num_threads <- ${THREADS}"
echo "pattern_f <- ${PATTERN}"
echo "trunc_f <- ${TRUNC_F}"
echo "trunc_r <- ${TRUNC_R}"

############################################################################################
#                                      Export output                                       #
############################################################################################

if [ $PRIMER_ATTACHED -eq 1 ]; then
  # remove files with trimmed sequences
  if [ -d "${INTERMEDIATE_DIR}/trimmed" ]; then
    rm -rf "${INTERMEDIATE_DIR}/trimmed"
  fi

  # remove files with untrimmed sequences
  if [ -d "${INTERMEDIATE_DIR}/untrimmed" ]; then
    rm -rf "${INTERMEDIATE_DIR}/untrimmed"
  fi
fi

# remove files with filtered sequences
if [ -d "${INTERMEDIATE_DIR}/filtered" ]; then
  rm -rf "${INTERMEDIATE_DIR}/filtered"
fi

# Reformat dna sequences and more to output folder
awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' "${INTERMEDIATE_DIR}/dna-sequences.fasta" > "${OUTPUT_DIR}/tmp.fasta"  && mv "${OUTPUT_DIR}/tmp.fasta" "${OUTPUT_DIR}/dna-sequences.fasta"

############################################################################################
#                                      Reformat output                                     #
############################################################################################

cp "${INTERMEDIATE_DIR}/dna-sequences.fasta" "${OUTPUT_DIR}/$(basename dna-sequences.fasta .fasta)_post.fasta"

awk 'BEGIN {FS=OFS="\t"} 
     NR==1 {for (i=1; i<=NF; i++) {if ($i != "ID") $i=tolower($i)}} 
     {print}' "${INTERMEDIATE_DIR}/taxonomy.tsv" > "${OUTPUT_DIR}/$(basename taxonomy.tsv .tsv)_post.tsv"

sed '1s/^/ID\t/' "${INTERMEDIATE_DIR}/feature-table.tsv" > "${OUTPUT_DIR}/$(basename feature-table.tsv .tsv)_post.tsv"

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
