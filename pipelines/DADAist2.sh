#!/bin/bash
#SBATCH -p rome
#SBATCH -t 2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/pipelines/DADAist2/log/DADAist2_output_%j.out

############################################################################################
#                              Pipeline configuration                                      #
############################################################################################

PIPELINE="DADAist2"

############################################################################################
#                                   Run configuration                                      #
############################################################################################

if [ -z "$SLURM_CPUS_PER_TASK" ]; then
  THREADS=96 
else
  THREADS=$SLURM_CPUS_PER_TASK
fi

# source conda environments
CONDA_BASE=$(conda info --base)
CONDA_SH="$CONDA_BASE/etc/profile.d/conda.sh"
source $CONDA_SH

# Set variables (which data to use, where to find it, and which primers are attached)
# Set default values
DEFAULT_DATASET="mock5"
DEFAULT_DATA_DIR="/gpfs/work2/0/gs/danique/mockrobiota/mock5/fastq_demux"
DEFAULT_FORWARD_PRIMER="GTGCCAGCMGCCGCGGTAA"
DEFAULT_REVERSE_PRIMER="GGACTACHVGGGTWTCTAAT"
DEFAULT_PRIMER_ATTACHED=0
DEFAULT_SEPARATOR="_"
DEFAULT_TRUNC_F=190
DEFAULT_TRUNC_R=135
DEFAULT_REMOVEOLD=0

# Set variables
while getopts ":d:i:f:r:a:s:t:u:y:h:" opt; do
  case $opt in
    d) DATASET="$OPTARG";;
    i) DATA_DIR="$OPTARG";;
    f) FORWARD_PRIMER="$OPTARG";;
    r) REVERSE_PRIMER="$OPTARG";;
    a) PRIMER_ATTACHED="$OPTARG";;
    s) SEPARATOR=="$OPTARG";;
    t) TRUNC_F="$OPTARG";;
    u) TRUNC_R="$OPTARG";; 
    y) REMOVEOLD="$OPTARG";;
    h)
      echo "Usage: \$0 [options]"
      echo "Options:"
      echo "  -d <dataset>            Specify dataset. This will be the directory name in which intermediate (intermediate/<DATASET>) and output (output/<DATASET>) files are stored"
      echo "  -i <data_dir>           Specify data directory. This directory should contain all fastq.gz files to be processed"
      echo "  -f <forward_primer>     Specify forward primer"
      echo "  -r <reverse_primer>     Specify reverse primer"
      echo "  -a <primer_attached>    Specify whether primers are still attached."
      echo "  -s <separator>          Specify which character is used to separate the sample ID from the rest of the fastq filename."         
      echo "  -t <truncate_forward> Specify at which position to truncate forward reads" 
      echo "  -u <truncate_reverse> Specify at which position to truncate reverse reads"
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
SEPARATOR=${SEPARATOR:-$DEFAULT_SEPARATOR}
TRUNC_F=${TRUNC_F:-$DEFAULT_TRUNC_F}
TRUNC_R=${TRUNC_R:-$DEFAULT_TRUNC_R}
REMOVEOLD=${REMOVEOLD:-$DEFAULT_REMOVEOLD}

# Pipeline Directories
BASE_DIR="/gpfs/work2/0/gs/danique/pipelines/pipelines/${PIPELINE}"
INTERMEDIATE_DIR="${BASE_DIR}/intermediate/${DATASET}"
OUTPUT_DIR="${BASE_DIR}/output/${DATASET}"
REF_DIR="${BASE_DIR}/intermediate"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${INTERMEDIATE_DIR}"

# Remove files from previous run if REMOVEOLD is set to 1
if [ $REMOVEOLD -eq 1 ]; then
  echo "removing intermediate files from previous run on ${DATASET}"
  rm -rf "${INTERMEDIATE_DIR}"/* # ensure dir is empty
fi

NEW_DATA_DIR="${INTERMEDIATE_DIR}/DATA" 
mkdir -p ${NEW_DATA_DIR}

# Add symbolic link to data folder
find "$DATA_DIR/" -maxdepth 1 -name "*.fastq.gz" -exec ln -sf {} "${NEW_DATA_DIR}" \;

############################################################################################
#                                      Run Pipeline                                        #
############################################################################################

conda activate dadaist2

# Check number of reads + read length per sample
seqfu stats -b -n ${NEW_DATA_DIR}/*_R1_*.fastq.gz --multiqc "${INTERMEDIATE_DIR}/${DATASET}_mqc.txt"

# Download reference database
dadaist2-getdb -d dada2-silva-138 -o "${REF_DIR}"

# Create metadata file
dadaist2-metadata -s ${SEPARATOR} -i "${NEW_DATA_DIR}" > "${BASE_DIR}/metadata.tsv"

head "${BASE_DIR}/metadata.tsv"

# Run pipeline
if [ "$PRIMER_ATTACHED" -eq 0 ]; then
  # If no primers attached
  dadaist2 \
    -i "${NEW_DATA_DIR}" \
    -o "${INTERMEDIATE_DIR}" \
    -m "${BASE_DIR}/metadata.tsv" \
    -d "${REF_DIR}/silva_nr_v138_train_set.fa.gz" \
    -no-trim \
    -s "${SEPARATOR}" \
    --trunc-len-1 "${TRUNC_F}" \
    --trunc-len-2 "${TRUNC_R}" \
    -t "${THREADS}"
else
  # If primers attached
  dadaist2 \
    -i "${NEW_DATA_DIR}" \
    -o "${INTERMEDIATE_DIR}" \
    -m "${BASE_DIR}/metadata.tsv" \
    -d "${REF_DIR}/silva_nr_v138_train_set.fa.gz" \
    -primers "${FORWARD_PRIMER}:${REVERSE_PRIMER}" \
    -s "${SEPARATOR}" \
    --trunc-len-1 "${TRUNC_F}" \
    --trunc-len-2 "${TRUNC_R}" \
    -t "${THREADS}"
fi

############################################################################################
#                                      Export output                                       #
############################################################################################

# Copy files to output dir
# Sequences
awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' "${INTERMEDIATE_DIR}/rep-seqs.fasta" > "${OUTPUT_DIR}/dna-sequences.fasta"

# Feature fable
cp "${INTERMEDIATE_DIR}/feature-table.tsv" "${OUTPUT_DIR}/feature-table.tsv"

# Taxonomy
cp "${INTERMEDIATE_DIR}/taxonomy.txt" "${OUTPUT_DIR}/taxonomy.tsv"

############################################################################################
#                                      Reformat output                                     #
############################################################################################

cp "${OUTPUT_DIR}/dna-sequences.fasta" "${OUTPUT_DIR}/$(basename dna-sequences.fasta .fasta)_post.fasta"
sed '1s/^#OTU //' "${OUTPUT_DIR}/feature-table.tsv" | \
awk 'BEGIN {FS=OFS="\t"} 
     NR==1 {for (i=1; i<=NF; i++) {split($i, a, "_"); $i=a[1]}} 
     {print}' > "${OUTPUT_DIR}/$(basename feature-table.tsv .tsv)_post.tsv"

sed 's/ /	/g' "${OUTPUT_DIR}/taxonomy.tsv" | \
sed '1s/^/ID\t/; 1!s/^/ASV/' | \
awk -F'\t' '
NR==1 {
    for (i=1; i<=NF; i++) header[$i] = i;
    print "ID", "kingdom", "phylum", "class", "order", "family", "genus";  # New header
    next
}
{
    print $header["ID"], $header["Kingdom"], $header["Phylum"], $header["Class"], $header["Order"], $header["Family"], $header["Genus"];
}' OFS='\t' > "${OUTPUT_DIR}/$(basename taxonomy.tsv .tsv)_post.tsv"

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