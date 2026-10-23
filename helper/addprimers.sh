#!/bin/bash
#SBATCH -p rome
#SBATCH -t 02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/log/sim/addprimers_%A_%a.out
#SBATCH --job-name=addprimers

# Adds the 515F/806R primer sequences, with maximum-quality scores, 
# to reads that don’t have them, for pipelines that require primers.

set -euo pipefail

# Expected input: addprimers.sh {DATA_DIR}
BASE_DIR=$1
cd "$BASE_DIR"

# New directory
mkdir -p "${BASE_DIR}/primer"

# Set primers
forward="GTGCCAGCAGCCGCGGTAA"
reverse="GGACTACAAGGGTATCTAAT"

# Calculate the length of the primer sequence
forward_length=${#forward}
reverse_length=${#reverse}

# Generate a string of maximum Phred quality scores (Phred+33, 'I' = 40)
forward_quality=$(printf 'I%.0s' $(seq 1 $forward_length))
reverse_quality=$(printf 'I%.0s' $(seq 1 $reverse_length))

# Export variables for GNU Parallel
export BASE_DIR forward reverse forward_quality reverse_quality

process_fastq() {
    ffile="$1"
    rfile="${ffile/R1/R2}"

    # Forward read
    pigz -dc "$ffile" | awk -v primer="$forward" -v qual="$forward_quality" '
        { if (NR % 4 == 1) print $0;
          else if (NR % 4 == 2) print primer $0;
          else if (NR % 4 == 3) print $0;
          else if (NR % 4 == 0) print qual $0;
        }' | pigz -c > "${BASE_DIR}/primer/$ffile"

    # Reverse read
    pigz -dc "$rfile" | awk -v primer="$reverse" -v qual="$reverse_quality" '
        { if (NR % 4 == 1) print $0;
          else if (NR % 4 == 2) print primer $0;
          else if (NR % 4 == 3) print $0;
          else if (NR % 4 == 0) print qual $0;
        }' | pigz -c > "${BASE_DIR}/primer/$rfile"
}
export -f process_fastq

# Find all R1 files and process in parallel (4 jobs at a time)
find . -maxdepth 1 -name '*R1*.fastq.gz' -print0 | \
    xargs -0 -n 1 -P 4 bash -c 'process_fastq "$0"' 

echo "All files processed"
