#!/bin/bash
#SBATCH -p rome
#SBATCH -t 5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/log/qc_output_%j.out
#SBATCH --mem=48G

if [ -z "$SLURM_CPUS_PER_TASK" ]; then
  THREADS=1 
else
  THREADS=$SLURM_CPUS_PER_TASK
fi

# source conda environments
CONDA_BASE=$(conda info --base)
CONDA_SH="$CONDA_BASE/etc/profile.d/conda.sh"
source $CONDA_SH

while getopts ":d:i:e:h:" opt; do
  case $opt in
    d) DATASET="$OPTARG";;
    i) DATA_DIR="$OPTARG";;
    e) EXT="$OPTARG";;
    h)
      echo "Usage: \$0 [options]"
      echo "Options:"
      echo "  -d <dataset>            Specify dataset name"
      echo "  -i <data_dir>           Specify data directory. This directory should contain all fastq.gz files to be processed"
      echo "  -e <extension>          Specify file extension for fastq files"
      echo "  -h                      Display this help message"
      exit 0
      ;;
    \?) echo "Invalid option: -$OPTARG"; exit 1;;
  esac
done

if [ -z "$EXT" ]; then
  EXT=".fastq.gz"
fi

OUTPUT_DIR="/gpfs/work2/0/gs/danique/pipelines/qc_output/${DATASET}"
mkdir -p $OUTPUT_DIR

# Empty output dir
rm -rf ${OUTPUT_DIR}/*

# activate conda environment
conda activate multiqc 

NUM_FILES=$(ls "$DATA_DIR"/*${EXT} | wc -l)
if [ "$NUM_FILES" -lt "$THREADS" ]; then
  THREADS=$NUM_FILES
fi

if [ "$NUM_FILES" -eq 0 ]; then
  echo "No ${EXT} files found in $DATA_DIR. Exiting."
  exit 1
else
  echo "Total number of ${EXT} files: $NUM_FILES"
fi

echo "Running FastQC on all ${EXT} files in $DATA_DIR..."
fastqc "$DATA_DIR"/*${EXT} -o "$OUTPUT_DIR" -t ${THREADS}

# Create forward.txt and reverse.txt files
if [ -f "${OUTPUT_DIR}/forward.txt" ]; then
  rm "${OUTPUT_DIR}/forward.txt"
fi

if [ -f "${OUTPUT_DIR}/reverse.txt" ]; then
  rm "${OUTPUT_DIR}/reverse.txt"
fi

touch ${OUTPUT_DIR}/forward.txt
touch ${OUTPUT_DIR}/reverse.txt

for file in "$OUTPUT_DIR"/*_R1_*fastqc.zip; do
  rfile=$(echo $file | sed 's/_R1_/_R2_/')
  echo $file >> ${OUTPUT_DIR}/forward.txt
  echo $rfile >> ${OUTPUT_DIR}/reverse.txt
done

# Run MultiQC to aggregate the FastQC results for reverse reads
echo "Running MultiQC to aggregate FastQC results for forward reads..."
multiqc -l ${OUTPUT_DIR}/forward.txt -n "$OUTPUT_DIR/forward" -f --clean-up
echo "Running MultiQC to aggregate FastQC results for reverse reads..."
multiqc -l ${OUTPUT_DIR}/reverse.txt -n "$OUTPUT_DIR/reverse" -f --clean-up

# Clean up
rm -rf ${OUTPUT_DIR}/*fastqc*