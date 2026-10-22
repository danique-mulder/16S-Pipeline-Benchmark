#!/bin/bash
#SBATCH -p rome
#SBATCH -t 2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/pipelines/nfcore_ampliseq/log/nfcore_ampliseq_output_%j.out

############################################################################################
#                              Pipeline configuration                                      #
############################################################################################

PIPELINE="nfcore_ampliseq"

############################################################################################
#                                   Run configuration                                      #
############################################################################################

if [ -z "$SLURM_CPUS_PER_TASK" ]; then
  THREADS=96
else
  # THREADS=$SLURM_CPUS_PER_TASK
  THREADS=$(($SLURM_NTASKS * $SLURM_CPUS_PER_TASK))
fi

if [ -z "$SLURM_MEM_PER_NODE" ]; then
  MEMORY="128GB"
else
  MEMORY=$OMP_MEMORY_PER_NODE
fi

# Set variables (which data to use, where to find it, and which primers are attached)
# Set default values
DEFAULT_DATASET="mock5"
DEFAULT_DATA_DIR="/gpfs/work2/0/gs/danique/mockrobiota/mock5/fastq_demux/"
DEFAULT_FORWARD_PRIMER="GTGCCAGCMGCCGCGGTAA"
DEFAULT_REVERSE_PRIMER="GGACTACHVGGGTWTCTAAT"
DEFAULT_PRIMER_ATTACHED=0
DEFAULT_SEPARATOR="_"
DEFAULT_REMOVEOLD=0

while getopts ":d:i:f:r:a:s:y:h:" opt; do
  case $opt in
    d) DATASET="$OPTARG";;
    i) DATA_DIR="$OPTARG";;
    f) FORWARD_PRIMER="$OPTARG";;
    r) REVERSE_PRIMER="$OPTARG";;
    a) PRIMER_ATTACHED="$OPTARG";;
    s) SEPARATOR=="$OPTARG";;
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
REMOVEOLD=${REMOVEOLD:-$DEFAULT_REMOVEOLD}

# Pipeline Directories
BASE_DIR="/gpfs/work2/0/gs/danique/pipelines/pipelines/${PIPELINE}"
INTERMEDIATE_DIR="${BASE_DIR}/intermediate/${DATASET}"
OUTPUT_DIR="${BASE_DIR}/output/${DATASET}"
NEXTFLOW="${BASE_DIR}/nextflow"

mkdir -p ${INTERMEDIATE_DIR}
mkdir -p ${OUTPUT_DIR}

# Remove files from previous run if REMOVEOLD is set to 1
if [ $REMOVEOLD -eq 1 ]; then
  echo "removing intermediate files from previous run on ${DATASET}"
  rm -rf "${INTERMEDIATE_DIR}"/* # ensure dir is empty
fi

############################################################################################
#                                   Create manifest_file                                   #
############################################################################################

echo -e "sampleID\tforwardReads\treverseReads" > "${INTERMEDIATE_DIR}/manifest.tsv" # File header

for filepath in ${DATA_DIR}/*_R1_*.fastq.gz; do 
  file=$(basename "$filepath")
  id=$(echo "$file" | awk -F"$SEPARATOR" '{print $1}')
  id="${id//./_}"
  path_forward=$filepath
  path_reverse="${path_forward//_R1_/_R2_}"
  
  echo -e "${id}\t${path_forward}\t${path_reverse}" >> "${INTERMEDIATE_DIR}/manifest.tsv"  
done 

head -n 3 "${INTERMEDIATE_DIR}/manifest.tsv" # Check manifest file

############################################################################################
#                                      Run Pipeline                                        #
############################################################################################
cd ${BASE_DIR}

# If primer, remove primers first (cutadapt)
if [ "$PRIMER_ATTACHED" -eq 1 ]; then
  $NEXTFLOW run nf-core/ampliseq \
     -profile singularity \
     -r 2.11.0 \
     --input "${INTERMEDIATE_DIR}/manifest.tsv" \
     --FW_primer $FORWARD_PRIMER \
     --RV_primer $REVERSE_PRIMER \
     --outdir $INTERMEDIATE_DIR \
     --max_cpus $THREADS \
     --max_memory $MEMORY
else
  $NEXTFLOW run nf-core/ampliseq \
       -profile singularity \
       -r 2.11.0 \
       --input "${INTERMEDIATE_DIR}/manifest.tsv" \
       --skip_cutadapt \
       --outdir $INTERMEDIATE_DIR \
       --max_cpus $THREADS \
       --max_memory $MEMORY
fi

############################################################################################
#                                      Export output                                       #
############################################################################################

awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' "${INTERMEDIATE_DIR}/dada2/ASV_seqs.fasta" > "${OUTPUT_DIR}/dna-sequences.fasta" # dna sequences
cp "${INTERMEDIATE_DIR}/dada2/ASV_table.tsv" "${OUTPUT_DIR}/feature-table.tsv" # feature table
cp "${INTERMEDIATE_DIR}/dada2/ASV_tax.silva_138.tsv" "${OUTPUT_DIR}/taxonomy.tsv" # taxonomy

############################################################################################
#                                      Reformat output                                     #
############################################################################################

# dna-sequences
cp "${OUTPUT_DIR}/dna-sequences.fasta" "${OUTPUT_DIR}/$(basename dna-sequences.fasta .fasta)_post.fasta"

# feature table
sed '1s/^ASV_//' "${OUTPUT_DIR}/feature-table.tsv" | \
awk 'BEGIN {FS=OFS="\t"} 
     NR==1 {for (i=1; i<=NF; i++) {split($i, a, "_"); $i=a[1]}} 
     {print}' > "${OUTPUT_DIR}/$(basename feature-table.tsv .tsv)_post.tsv" 

# taxonomy
cut -f1-$(($(awk -F'\t' '{print NF; exit}' "${OUTPUT_DIR}/taxonomy.tsv") - 3)) "${OUTPUT_DIR}/taxonomy.tsv" | \
awk 'BEGIN {OFS="\t"; print "ID", "kingdom", "phylum", "class", "order", "family", "genus", "species"} {if (NR > 1) print}' > "${OUTPUT_DIR}/$(basename taxonomy.tsv .tsv)_post.tsv" 

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
