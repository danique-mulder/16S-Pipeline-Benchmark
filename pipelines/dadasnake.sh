#!/bin/bash
#SBATCH -p rome
#SBATCH -t 2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/pipelines/dadasnake/log/dadasnake_output_%j.out

############################################################################################
#                              Pipeline configuration                                      #
############################################################################################

PIPELINE="dadasnake"
RDP="https://mothur.s3.us-east-2.amazonaws.com/wiki/trainset9_032012.pds.zip"

############################################################################################
#                                   Run configuration                                      #
############################################################################################

if [ -z "$SLURM_CPUS_PER_TASK" ]; then
  THREADS=96
else
  # THREADS=$SLURM_CPUS_PER_TASK
  THREADS=$(($SLURM_NTASKS * $SLURM_CPUS_PER_TASK))
fi

# Set variables (which data to use, where to find it, and which primers are attached)
# Set default values
DEFAULT_DATASET="mock5"
DEFAULT_DATA_DIR="/gpfs/work2/0/gs/danique/mockrobiota/mock5/fastq_demux"
DEFAULT_FORWARD_PRIMER="GTGCCAGCMGCCGCGGTAA"
DEFAULT_REVERSE_PRIMER="GGACTACHVGGGTWTCTAAT"
DEFAULT_PRIMER_ATTACHED=0
DEFAULT_SEPARATOR="_"
DEFAULT_MACHINE="miSeq"
DEFAULT_TRUNC_F=190
DEFAULT_TRUNC_R=135
DEFAULT_REMOVEOLD=0

while getopts ":d:i:f:r:a:s:m:t:u:y:h:" opt; do
  case $opt in
    d) DATASET="$OPTARG";;
    i) DATA_DIR="$OPTARG";;
    f) FORWARD_PRIMER="$OPTARG";;
    r) REVERSE_PRIMER="$OPTARG";;
    a) PRIMER_ATTACHED="$OPTARG";;
    s) SEPARATOR=="$OPTARG";;
    m) MACHINE=="$OPTARG";;
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
      echo "  -m <machine>            Specify which machine generated the sequences (NovaSeq/MiSeq/HiSeq.)"
      echo "  -t <truncate_forward>   Specify at which position to truncate forward reads" 
      echo "  -u <truncate_reverse>   Specify at which position to truncate reverse reads" 
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
MACHINE=${MACHINE:-$DEFAULT_MACHINE}
TRUNC_F=${TRUNC_F:-$DEFAULT_TRUNC_F}
TRUNC_R=${TRUNC_R:-$DEFAULT_TRUNC_R}
REMOVEOLD=${REMOVEOLD:-$DEFAULT_REMOVEOLD}

# Pipeline Directories
BASE_DIR="/gpfs/work2/0/gs/danique/pipelines/pipelines/${PIPELINE}"
INTERMEDIATE_DIR="${BASE_DIR}/intermediate/${DATASET}"
OUTPUT_DIR="${BASE_DIR}/output/${DATASET}"
DADASNAKE_DIR="${BASE_DIR}/dadasnake"
DADASNAKE="${BASE_DIR}/dadasnake/dadasnake"
REF_DIR="${BASE_DIR}/intermediate/"

mkdir -p $INTERMEDIATE_DIR
mkdir -p $OUTPUT_DIR

# Remove files from previous run if REMOVEOLD is set to 1
if [ $REMOVEOLD -eq 1 ]; then
  echo "removing intermediate files from previous run on ${DATASET}"
  rm -rf "${INTERMEDIATE_DIR}"/* # ensure dir is empty
fi

############################################################################################
#                                   Install databases                                      #
############################################################################################

if [ ! -f "${REF_DIR}/trainset9_032012.pds.taxonomy" ]; then
  wget -P "$REF_DIR/" $RDP
  unzip "${REF_DIR}/trainset9_032012.pds.zip" -d $REF_DIR
  mv "${REF_DIR}/trainset9_032012.pds.tax" "${REF_DIR}/trainset9_032012.pds.taxonomy" 
fi

############################################################################################
#                                Create manifest_file                                      #
############################################################################################

# File header
echo -e "library\tr1_file\tr2_file" > "${INTERMEDIATE_DIR}/manifest.tsv"

for filepath in ${DATA_DIR}/*_R1_*.fastq.gz; do 

  file=$(basename "$filepath")
  id=$(echo "$file" | awk -F"$SEPARATOR" '{print $1}')
  id="${id//./_}"
  
  path_forward=$file
  path_reverse="${path_forward//_R1_/_R2_}"

  echo -e "${id}\t${path_forward}\t${path_reverse}" >> "${INTERMEDIATE_DIR}/manifest.tsv" 
done 

head -n 3 "${INTERMEDIATE_DIR}/manifest.tsv"

############################################################################################
#                                 Adjust configuration file                                #
############################################################################################

# Specify dadasnake variable configuration file
cat <<EOL > "${DADASNAKE_DIR}/VARIABLE_CONFIG.sbatch"
SNAKEMAKE_VIA_CONDA,true
SNAKEMAKE_EXTRA_ARGUMENTS,
LOADING_MODULES,
SUBMIT_COMMAND,sbatch
BIND_JOBS_TO_MAIN,false
NODENAME_VAR,
SCHEDULER,slurm
MAX_THREADS,20
BIGMEM_CORES,10	
BIGMEM_MEM_EACH,6G
NORMAL_MEM_EACH,6G
LOCK_SETTINGS,false		
EOL


sed 's/,/\t/g' "${DADASNAKE_DIR}/VARIABLE_CONFIG.sbatch" > "${DADASNAKE_DIR}/VARIABLE_CONFIG"
cp "${DADASNAKE_DIR}/auxiliary_files/dadasnake_allSubmit" dadasnake
chmod 755 dadasnake

# Pipeline run configuration (yaml)
CONFIG_FILE="${DADASNAKE_DIR}/config/config.16S.${DATASET}.yaml"
cp "${DADASNAKE_DIR}/config/config.16S.yaml" ${CONFIG_FILE}

# Data specification
yq --arg DATA_DIR "$DATA_DIR" '.raw_directory = $DATA_DIR' -i -y ${CONFIG_FILE}
yq --arg INTERMEDIATE_DIR "$INTERMEDIATE_DIR" '.sample_table = $INTERMEDIATE_DIR + "/manifest.tsv"' -i -y ${CONFIG_FILE}
yq --arg INTERMEDIATE_DIR "$INTERMEDIATE_DIR" '.outputdir = $INTERMEDIATE_DIR' -i -y ${CONFIG_FILE}
yq --arg INTERMEDIATE_DIR "$INTERMEDIATE_DIR" '.tmp_dir = $INTERMEDIATE_DIR' -i -y ${CONFIG_FILE}

# Primers
if [ "$PRIMER_ATTACHED" -eq 0 ]; then
  yq '.do_primers = false' -i -y ${CONFIG_FILE}
else
  yq --arg FORWARD_PRIMER "$FORWARD_PRIMER" '.primers.fwd.sequence = $FORWARD_PRIMER' -i -y ${CONFIG_FILE}
  yq --arg REVERSE_PRIMER "$REVERSE_PRIMER" '.primers.rvs.sequence = $REVERSE_PRIMER' -i -y ${CONFIG_FILE}
fi

# MACHINE
if [ "$MACHINE" = "novaSeq" ]; then
  yq '.nextseq_novaseq = true' -i -y ${CONFIG_FILE}
fi

# Trimming (change default to be suitable for input data)
yq --arg TRUNC_F "$TRUNC_F" '.filtering.trunc_length.fwd = ($TRUNC_F | tonumber)' -i -y ${CONFIG_FILE}
yq --arg TRUNC_R "$TRUNC_R" '.filtering.trunc_length.rvs = ($TRUNC_R | tonumber)' -i -y ${CONFIG_FILE}

# Taxonomic assignment
yq --arg REF_DIR "$REF_DIR" '.taxonomy.mothur.db_path = $REF_DIR' -i -y ${CONFIG_FILE}
yq '.taxonomy.mothur.tax_db = "trainset9_032012.pds"' -i -y ${CONFIG_FILE}
yq '.taxonomy.mothur.db_short_names = "trainset9_032012.pds"' -i -y ${CONFIG_FILE}

# Blast and postprocessing
yq '.blast.do = false' -i -y ${CONFIG_FILE}
yq '.postprocessing.picrust2.do = false' -i -y ${CONFIG_FILE}

############################################################################################
#                                      Run Pipeline                                        #
############################################################################################
cd $DADASNAKE_DIR
$DADASNAKE -l -r -u -n ${DATASET} ${CONFIG_FILE}
$DADASNAKE -t ${THREADS} -l -r -n ${DATASET} ${CONFIG_FILE} 

############################################################################################
#                                      Export output                                       #
############################################################################################

# DNA sequences
awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' "${INTERMEDIATE_DIR}/sequenceTables/all.seqs.fasta" > "${OUTPUT_DIR}/dna-sequences.fasta"

# Feature table
cp "${INTERMEDIATE_DIR}/sequenceTables/all.seqTab.tsv" "${OUTPUT_DIR}/feature-table.tsv"

# Taxonomy
cp "${INTERMEDIATE_DIR}/sequenceTables/tax.mothur.trainset9_032012.pds.tsv" "${OUTPUT_DIR}/taxonomy.tsv" 

# With 97% post-ASV clustering
# DNA sequences, but print the whole sequence on one line (otu)
awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' "${INTERMEDIATE_DIR}/clusteredTables/consensus.fasta" > "${OUTPUT_DIR}/dna-sequences_otu.fasta"

# Feature table (otu)
cp "${INTERMEDIATE_DIR}/clusteredTables/clusteredTab.tsv" "${OUTPUT_DIR}/feature-table_otu.tsv"

# Taxonomy (otu)
cp "${INTERMEDIATE_DIR}/clusteredTables/tax.mothur.trainset9_032012.pds.tsv" "${OUTPUT_DIR}/taxonomy_otu.tsv" 

############################################################################################
#                                      Reformat output                                     #
############################################################################################

# dna-sequences
cp "${OUTPUT_DIR}/dna-sequences.fasta" "${OUTPUT_DIR}/$(basename dna-sequences.fasta .fasta)_post.fasta"
cp "${OUTPUT_DIR}/dna-sequences_otu.fasta" "${OUTPUT_DIR}/$(basename dna-sequences_otu.fasta .fasta)_post.fasta"

# feature table (asv)
cut -f2-  "${OUTPUT_DIR}/feature-table.tsv" | \
awk 'BEGIN {OFS=FS="\t"} {print $NF, $0}' | \
awk -F'\t' '{ $NF=""; print $0 }' OFS='\t' | \
sed '1s/^ASV/ID/' > "${OUTPUT_DIR}/$(basename feature-table.tsv .tsv)_post.tsv"

# feature table (otu)
cut -f1-$(($(awk -F'\t' '{print NF; exit}' "${OUTPUT_DIR}/feature-table_otu.tsv") - 10)) "${OUTPUT_DIR}/feature-table_otu.tsv" | \
cut -f2- | cut -f1,3- | \
sed '1s/^OTU/ID/' > "${OUTPUT_DIR}/$(basename feature-table_otu.tsv .tsv)_post.tsv"

# Taxonomy (asv)
sed 's/;/\t/g;s/"//g;s/([^(]*)//g' "${OUTPUT_DIR}/taxonomy.tsv" | \
awk 'BEGIN {OFS="\t"; print "ID", "kingdom", "phylum", "class", "order", "family", "genus"} {print}' > "${OUTPUT_DIR}/$(basename taxonomy.tsv .tsv)_post.tsv"

sed 's/;/\t/g;s/"//g;s/([^(]*)//g' "${OUTPUT_DIR}/taxonomy_otu.tsv" | \
awk 'BEGIN {OFS="\t"; print "ID", "kingdom", "phylum", "class", "order", "family", "genus"} {print}' > "${OUTPUT_DIR}/$(basename taxonomy_otu.tsv .tsv)_post.tsv"

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
numgenus=$(awk -F'\t' '{print $7 }' "${OUTPUT_DIR}/taxonomy_otu_post.tsv" | sort | uniq -c | wc -l)

echo "asv: ${numasv}" > "${OUTPUT_DIR}/stat_otu.txt"
echo "genera: ${numgenus}" >> "${OUTPUT_DIR}/stat_otu.txt"

awk -F'\t' '{print $7 }' "${OUTPUT_DIR}/taxonomy_otu_post.tsv" | sort | uniq -c | while read line; do
  echo "$line" >> "${OUTPUT_DIR}/stat_otu.txt"
done