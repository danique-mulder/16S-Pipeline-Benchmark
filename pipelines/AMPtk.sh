#!/bin/bash
#SBATCH -p rome
#SBATCH -t 2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/pipelines/AMPtk/log/AMPtk_output_%j.out

############################################################################################
#                              Pipeline configuration                                      #
############################################################################################

PIPELINE="AMPtk"

############################################################################################
#                                   Run configuration                                      #
############################################################################################

if [ -z "$SLURM_CPUS_PER_TASK" ]; then
  THREADS=96 
else
  # THREADS=$SLURM_CPUS_PER_TASK
  THREADS=$(($SLURM_NTASKS * $SLURM_CPUS_PER_TASK))
fi

# Access r libraries
R_LIBS=$HOME/rpackages:$R_LIBS
export R_LIBS

# Load R module (pre-installed in Snellius)
module load 2023
module load R/4.3.2-gfbf-2023a

# source conda environments
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
      echo "  -y <remove_old>         Specify whether to remove old files (0 = no, 1 = yes)"      
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

# Pipeline Directories
BASE_DIR="/gpfs/work2/0/gs/danique/pipelines/pipelines/${PIPELINE}"
INTERMEDIATE_DIR="${BASE_DIR}/intermediate/${DATASET}"
OUTPUT_DIR="${BASE_DIR}/output/${DATASET}"

SCRIPT_DIR="/gpfs/work2/0/gs/danique/pipelines/scripts"
LOG_DIR="${BASE_DIR}/log"
REF_DIR="${BASE_DIR}/intermediate/"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${INTERMEDIATE_DIR}"

# Remove files from previous run if REMOVEOLD is set to 1
if [ $REMOVEOLD -eq 1 ]; then
  echo "removing intermediate files from previous run on ${DATASET}"
  rm -rf "${INTERMEDIATE_DIR}"/* # ensure dir is empty
fi

############################################################################################
#                                      Run Pipeline                                        #
############################################################################################

cd $INTERMEDIATE_DIR

# Activate conda environment with the installation
conda activate amptk

# Run pre-processing, either removing the primer if still attached or not if not attached
if [ "$PRIMER_ATTACHED" -eq 1 ]; then
  amptk illumina -i "${DATA_DIR}" -o "${DATASET}" -f "$FORWARD_PRIMER" -r "$REVERSE_PRIMER" -l 275  --cpus 8 --merge_method vsearch --cleanup
else
  amptk illumina -i "${DATA_DIR}" -o "${DATASET}" -f "" -r "" -l 275 --cpus 8 --merge_method vsearch --cleanup
fi

# Clustering data > there is no default for truncation (-l argument; truncation after merging)
amptk dada2 -i "${DATASET}.demux.fq.gz" -o "${DATASET}" --platform illumina -l 250

# Filter to detect index-bleeding
amptk filter -i "${DATASET}.otu_table.txt" -f "${DATASET}.ASVs.fa" -o "${DATASET}" --show_stats 
amptk filter -i "${DATASET}.cluster.otu_table.txt" -f "${DATASET}.cluster.otus.fa" -o "${DATASET}.cluster" --show_stats 

# LULU2 post-clustering
# LULU2 cannot be ran on a single sample. Check if there are multiple samples. If so, run LULU2
if [ $(ls "${DATA_DIR}"/*R1*.fastq.gz | wc -l) -gt 1 ]; then
  amptk lulu -i "${DATASET}.final.txt" -f "${DATASET}.filtered.otus.fa" -o "${DATASET}"
  amptk lulu -i "${DATASET}.cluster.final.txt" -f "${DATASET}.cluster.filtered.otus.fa" -o "${DATASET}.cluster"

  amptk install -i 16S
  amptk taxonomy -f "${DATASET}.lulu.otus.fa" -i "${DATASET}.lulu.otu_table.txt" -o "${DATASET}" -m "${DATASET}.mapping_file.txt" -d 16S --method hybrid
  amptk taxonomy -f "${DATASET}.cluster.lulu.otus.fa" -i "${DATASET}.cluster.lulu.otu_table.txt" -o "${DATASET}.cluster" -m "${DATASET}.mapping_file.txt" -d 16S --method hybrid
else 
  amptk install -i 16S
  amptk taxonomy -f "${DATASET}.ASVs.fa" -i "${DATASET}.otu_table.txt" -o "${DATASET}" -m "${DATASET}.mapping_file.txt" -d 16S --method hybrid
  amptk taxonomy -f "${DATASET}.cluster.otus.fa" -i "${DATASET}.cluster.otu_table.txt" -o "${DATASET}.cluster" -m "${DATASET}.mapping_file.txt" -d 16S --method hybrid
fi

conda deactivate

############################################################################################
#                                      Export output                                       #
############################################################################################

# Clean up directory
mv ${INTERMEDIATE_DIR}/*.log $LOG_DIR

# If LULU2 was ran, use the LULU2 output files, otherwise use the regular output files
if [ $(ls "${DATA_DIR}"/*R1*.fastq.gz | wc -l) -gt 1 ]; then
  # Feature table
  cp "${DATASET}.lulu.otu_table.txt" "${OUTPUT_DIR}/feature-table.tsv"
  cp "${DATASET}.cluster.lulu.otu_table.txt" "${OUTPUT_DIR}/feature-table_otu.tsv"
else
  # Feature table
  cp "${DATASET}.otu_table.txt" "${OUTPUT_DIR}/feature-table.tsv"
  cp "${DATASET}.cluster.otu_table.txt" "${OUTPUT_DIR}/feature-table_otu.tsv"
fi

# Move taxonomy, dna sequences and otutable to output file
# DNA sequences, but print the whole sequence on one line
awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' "${DATASET}.cluster.otus.taxonomy.fa" > "${OUTPUT_DIR}/dna-sequences_otu.fasta"
awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' "${DATASET}.otus.taxonomy.fa" > "${OUTPUT_DIR}/dna-sequences.fasta"

# Taxonomic assignment
cp "${DATASET}.taxonomy.txt" "${OUTPUT_DIR}/taxonomy.tsv"
cp "${DATASET}.cluster.taxonomy.txt" "${OUTPUT_DIR}/taxonomy_otu.tsv"

############################################################################################
#                                      Reformat output                                     #
############################################################################################

# dna-sequences (asv)
awk '{if($0 ~ /^>/) print $1; else print $0}' "${OUTPUT_DIR}/dna-sequences.fasta" > "${OUTPUT_DIR}/$(basename dna-sequences.fasta .fasta)_post.fasta"

# dna-sequences (otu)
awk '{if($0 ~ /^>/) print $1; else print $0}' "${OUTPUT_DIR}/dna-sequences_otu.fasta" > "${OUTPUT_DIR}/$(basename dna-sequences_otu.fasta .fasta)_post.fasta"

# feature table
sed '1s/#OTU //g' "${OUTPUT_DIR}/feature-table.tsv" > "${OUTPUT_DIR}/$(basename feature-table.tsv .tsv)_post.tsv" # asv
sed '1s/#OTU //g' "${OUTPUT_DIR}/feature-table_otu.tsv" > "${OUTPUT_DIR}/$(basename feature-table_otu.tsv .tsv)_post.tsv" #otu

# taxonomy (av
# taxonomy (otu)
cut -f1-$(($(awk -F'\t' '{print NF; exit}' "${OUTPUT_DIR}/taxonomy_otu.tsv") - 3)) "${OUTPUT_DIR}/taxonomy_otu.tsv" | \
sed 's/;/\t/g;s/,/\t/g;s/.\?://g;s/"//g' | cut -f1,3- | \
awk 'BEGIN {OFS="\t"; print "ID", "kingdom", "phylum", "class", "order", "family", "genus"} {if (NR > 1) print}' > "${OUTPUT_DIR}/$(basename taxonomy_otu.tsv .tsv)_post.tsv"

############################################################################################
#                                      Descriptive stats                                   #
############################################################################################

## ASV
# Write basic descriptive stats to txt file
numasv=$(awk '/^>/' "${OUTPUT_DIR}/dna-sequences_post.fasta" | wc -l) 
numgenus=$(awk -F'\t' '{print $7 }' "${OUTPUT_DIR}/taxonomy_post.tsv" | sort | uniq -c | wc -l)

echo "asv: ${numasv}" > "${OUTPUT_DIR}/stat.txt"
echo "genera: ${numgenus}" >> "${OUTPUT_DIR}/stat.txt"

awk -F'\t' '{print $7 }' "${OUTPUT_DIR}/taxonomy_post.tsv" | sort | uniq -c | while read line; do
  echo "$line" >> "${OUTPUT_DIR}/stat.txt"
done

## OTU (LULU)
# Write basic descriptive stats to txt file
numasv=$(awk '/^>/' "${OUTPUT_DIR}/dna-sequences_otu_post.fasta" | wc -l) 
numgenus=$(awk -F'\t' '{print $7 }' "${OUTPUT_DIR}/taxonomy_post.tsv" | sort | uniq -c | wc -l)

echo "asv: ${numasv}" > "${OUTPUT_DIR}/stat_otu.txt"
echo "genera: ${numgenus}" >> "${OUTPUT_DIR}/stat_otu.txt"

awk -F'\t' '{print $7 }' "${OUTPUT_DIR}/taxonomy_otu_post.tsv" | sort | uniq -c | while read line; do
  echo "$line" >> "${OUTPUT_DIR}/stat_otu.txt"
done