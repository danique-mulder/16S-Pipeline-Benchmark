#!/bin/bash
#SBATCH -p rome
#SBATCH -t 20:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/pipelines/LotuS2/log/LotuS2_output_%j.out

############################################################################################
#                              Pipeline configuration                                      #
############################################################################################

PIPELINE="LotuS2"

############################################################################################
#                                   Run configuration                                      #
############################################################################################

if [ -z "$SLURM_CPUS_PER_TASK" ]; then
  THREADS=96  # 128 GB in MB
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

# Set variables (which data to use, where to find it, and which primers are attached)
# Set default values
DEFAULT_DATASET="mock5" 
DEFAULT_DATA_DIR="/gpfs/work2/0/gs/danique/mockrobiota/mock5/fastq_demux/primer" 
DEFAULT_FORWARD_PRIMER="GTGCCAGCMGCCGCGGTAA" 
DEFAULT_REVERSE_PRIMER="GGACTACHVGGGTWTCTAAT" 
DEFAULT_PRIMER_ATTACHED=1 
DEFAULT_MACHINE="miSeq" 
DEFAULT_SEQUENCING_RUN=0 
DEFAULT_CLUSTER="asv"

while getopts ":d:i:f:r:a:m:l:c:h:" opt; do
  case $opt in
    d) DATASET="$OPTARG";;
    i) DATA_DIR="$OPTARG";;
    f) FORWARD_PRIMER="$OPTARG";;
    r) REVERSE_PRIMER="$OPTARG";;
    a) PRIMER_ATTACHED="$OPTARG";;
    m) MACHINE="$OPTARG";;
    l) SEQUENCING_RUN="$OPTARG";;
    c) CLUSTER="$OPTARG";;
    h)
      echo "Usage: \$0 [options]"
      echo "Options:"
      echo "  -d <dataset>         Specify dataset. This will be the directory name in which intermediate (intermediate/<DATASET>) and output (output/<DATASET>) files are stored"
      echo "  -i <data_dir>        Specify data directory. This directory should contain all fastq.gz files to be processed"
      echo "  -f <forward_primer>  Specify forward primer"
      echo "  -r <reverse_primer>  Specify reverse primer"
      echo "  -a <primer_attached> Specify whether primers are still attached."
      echo "  -m <machine>         Specify sequencing machine (miSeq/novaSeq). This determines the SDM configuration file."
      echo "  -l <sequencing_run>  Specify whether sequences were processed in different sequencing runs (0/1). If 1, DADA2 is ran separately for each run"
      echo "  -c <clustering>      Specify which clustering option to use: 1) asv = dada2 or 2) otu = uparse"     
      echo "  -h                   Display this help message"
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
MACHINE=${MACHINE:-$DEFAULT_MACHINE}
SEQUENCING_RUN=${SEQUENCING_RUN:-$DEFAULT_SEQUENCING_RUN}
CLUSTER=${CLUSTER:-$DEFAULT_CLUSTER}

# Pipeline Directories
BASE_DIR="/gpfs/work2/0/gs/danique/pipelines/pipelines/${PIPELINE}"
INTERMEDIATE_DIR="${BASE_DIR}/intermediate/${DATASET}/${CLUSTER}"
OUTPUT_DIR="${BASE_DIR}/output/${DATASET}"

mkdir -p "${INTERMEDIATE_DIR}"
mkdir -p "${OUTPUT_DIR}"

#rm -rf "${INTERMEDIATE_DIR}/"* # Lotus automatically overwrites previous files

LOTUS_DIR="${BASE_DIR}/lotus2"
LOTUS="${LOTUS_DIR}/lotus2" # lotus2 executable

############################################################################################
#                                       Create mapping file                                #
############################################################################################

# Move to lotus executable directory
cd "$LOTUS_DIR"

# Create mapping file
MAPPING_FILE="${INTERMEDIATE_DIR}/manifest.tsv"

if [ $PRIMER_ATTACHED -eq 0 ]; then
  echo "ERROR - Sequences need to have primers attached"
  exit 1
fi

# mapping file MUST start with #
echo -e "#SampleID\tfastqFile" > $MAPPING_FILE # File header

for filepath in ${DATA_DIR}/*_R1_*.fastq.gz; do 
  file=$(basename "$filepath")
  id=$(echo "$file" | awk -F"_" '{print $1}')
  path_forward=$file
  path_reverse="${path_forward//_R1_/_R2_}"
  
  path_total=$(echo "${path_forward},${path_reverse}")
  
  echo -e "${id}\t${path_total}" >> $MAPPING_FILE  
done 

head -n 10 $MAPPING_FILE # preview the mapping file

############################################################################################
#                                 Adjust configuration file                                #
############################################################################################

# Print LotuS2 Version number
$LOTUS -v

SDM="${INTERMEDIATE_DIR}/config.txt"
cp "${LOTUS_DIR}/configs/sdm_miSeq.txt" ${SDM} # same for NovaSeq

sed -i 's/maxBarcodeErrs\t0/maxBarcodeErrs\t-1/' "${SDM}" # No barcodes present

# Dynamically change minimum sequence length to length - 50 for shorter amplicons
seq=$(zcat "$(ls -1 "${DATA_DIR}"/*.fastq.gz | head -n 1)" | head -n 2 | tail -n 1)
seqlength=${#seq} 
minseqlength=$((seqlength - 50)) # Set minimum sequence length to 50
if [ $minseqlength -lt 200 ]; then
  sed -i "s/minSeqLength\t200/minSeqLength\t${minseqlength}/" "${SDM}" # Set minimum sequence length to the calculated value
fi

# sed -i 's/maxPrimerErrs\t0/maxPrimerErrs\t1/' "${SDM}" # Allow 1 primer error
# sed -i 's/PEheaderPairFmt  1/PEheaderPairFmt  0/' "${SDM}" # Does not run when set to default (1)


############################################################################################
#                                      Run Pipeline                                        #
############################################################################################

if [ $CLUSTER == "asv" ]; then 
  # DADA2
  $LOTUS \
  -i "${DATA_DIR}" \
  -o "${INTERMEDIATE_DIR}/results" \
  -m "${MAPPING_FILE}" \
  -forwardPrimer $FORWARD_PRIMER \
  -reversePrimer $REVERSE_PRIMER \
  -verbosity 3 \
  -CL dada2 \
  -ITSx 0 \
  -buildPhylo 0 \
  -lulu 1 \
  -p $MACHINE \
  -s $SDM \
  -t $THREADS
elif [ $CLUSTER == "otu" ]; then
  # UPARSE
  $LOTUS \
  -i "${DATA_DIR}" \
  -o "${INTERMEDIATE_DIR}/results" \
  -m "${MAPPING_FILE}" \
  -forwardPrimer $FORWARD_PRIMER \
  -reversePrimer $REVERSE_PRIMER \
  -verbosity 3 \
  -CL uparse \
  -ITSx 0 \
  -buildPhylo 0 \
  -lulu 1 \
  -p $MACHINE \
  -s $SDM \
  -t $THREADS
else
  echo "-c should be either 'asv' or 'otu'"
  exit 1
fi

############################################################################################
#                                      Export output                                       #
############################################################################################
 
# Copy files to the output directory
if [ $CLUSTER == "asv" ]; then  
  
  awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' "${INTERMEDIATE_DIR}/results/OTU.fna" > "${OUTPUT_DIR}/dna-sequences.fasta" 
  cp "${INTERMEDIATE_DIR}/results/OTU.txt" "${OUTPUT_DIR}/feature-table.tsv"
  cp "${INTERMEDIATE_DIR}/results/hiera_RDP.txt" "${OUTPUT_DIR}/taxonomy.tsv"
  
elif [ $CLUSTER == "otu" ]; then
  
  awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' "${INTERMEDIATE_DIR}/results/OTU.fna" > "${OUTPUT_DIR}/dna-sequences_otu.fasta" 
  cp "${INTERMEDIATE_DIR}/results/OTU.txt" "${OUTPUT_DIR}/feature-table_otu.tsv"
  cp "${INTERMEDIATE_DIR}/results/hiera_RDP.txt" "${OUTPUT_DIR}/taxonomy_otu.tsv"
  
fi

############################################################################################
#                                      Reformat output                                     #
############################################################################################

if [ $CLUSTER == "asv" ]; then 
  # dna-sequences
  cp "${OUTPUT_DIR}/dna-sequences.fasta" "${OUTPUT_DIR}/$(basename dna-sequences.fasta .fasta)_post.fasta" # asv
  
  # feature table
  awk 'NR==1 {sub(/^\t|^OTU\t/, "ID\t"); print; next} {print}' "${OUTPUT_DIR}/feature-table.tsv" > "${OUTPUT_DIR}/$(basename feature-table.tsv .tsv)_post.tsv"  # asv
  
  # taxonomy
  awk 'BEGIN {OFS="\t"; print "ID", "kingdom", "phylum", "class", "order", "family", "genus", "species"} {if (NR > 1) print}' "${OUTPUT_DIR}/taxonomy.tsv" > "${OUTPUT_DIR}/$(basename taxonomy.tsv .tsv)_post.tsv"  # asv

elif [ $CLUSTER == "otu" ]; then
  
  # dna-sequences
  cp "${OUTPUT_DIR}/dna-sequences_otu.fasta" "${OUTPUT_DIR}/$(basename dna-sequences_otu.fasta .fasta)_post.fasta" # otu
  
  # featue table
  awk 'NR==1 {sub(/^\t/, "ID\t"); print; next} {print}' "${OUTPUT_DIR}/feature-table_otu.tsv" > "${OUTPUT_DIR}/$(basename feature-table_otu.tsv .tsv)_post.tsv" # otu
  
  # taxonomy
  awk 'BEGIN {OFS="\t"; print "ID", "kingdom", "phylum", "class", "order", "family", "genus", "species"} {if (NR > 1) print}' "${OUTPUT_DIR}/taxonomy_otu.tsv" > "${OUTPUT_DIR}/$(basename taxonomy_otu.tsv .tsv)_post.tsv" # otu
  
fi

############################################################################################
#                                      Descriptive stats                                   #
############################################################################################

if [ $CLUSTER == "asv" ]; then  
  ## ASV
  # Write basic descriptive stats to txt file
  numasv=$(awk '/^>/' "${OUTPUT_DIR}/dna-sequences_post.fasta" | wc -l) 
  numgenus=$(awk -F'\t' '{print $7 }' "${OUTPUT_DIR}/taxonomy_post.tsv" | sort | uniq -c | wc -l)

  echo "asv: ${numasv}" > "${OUTPUT_DIR}/stat.txt"
  echo "genera: ${numgenus}" >> "${OUTPUT_DIR}/stat.txt"

  awk -F'\t' '{print $7 }' "${OUTPUT_DIR}/taxonomy_post.tsv" | sort | uniq -c | while read line; do
    echo "$line" >> "${OUTPUT_DIR}/stat.txt"
  done
elif [ $CLUSTER == "otu" ]; then
  ## OTU (LULU)
  # Write basic descriptive stats to txt file
  numasv=$(awk '/^>/' "${OUTPUT_DIR}/dna-sequences_otu_post.fasta" | wc -l) 
  numgenus=$(awk -F'\t' '{print $7 }' "${OUTPUT_DIR}/taxonomy_otu_post.tsv" | sort | uniq -c | wc -l)

  echo "asv: ${numasv}" > "${OUTPUT_DIR}/stat_otu.txt"
  echo "genera: ${numgenus}" >> "${OUTPUT_DIR}/stat_otu.txt"

  awk -F'\t' '{print $7 }' "${OUTPUT_DIR}/taxonomy_otu_post.tsv" | sort | uniq -c | while read line; do
    echo "$line" >> "${OUTPUT_DIR}/stat_otu.txt"
  done
fi
