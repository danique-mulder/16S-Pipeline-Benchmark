#!/bin/bash
#SBATCH -p rome
#SBATCH -t 2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/pipelines/MetaWorks/log/MetaWorks_output_%j.out

############################################################################################
#                              Pipeline configuration                                      #
############################################################################################

PIPELINE="MetaWorks"

############################################################################################
#                                   Run configuration                                      #
############################################################################################

if [ -z "$SLURM_CPUS_PER_TASK" ]; then
  THREADS=96  # 128 GB in MB
else
  # THREADS=$SLURM_CPUS_PER_TASK
  THREADS=$(($SLURM_NTASKS * $SLURM_CPUS_PER_TASK))
fi

# source conda environments
CONDA_BASE=$(conda info --base)
CONDA_SH="$CONDA_BASE/etc/profile.d/conda.sh"
source $CONDA_SH

# Set variables (which data to use, where to find it, and which primers are attached)
# Set default values
DEFAULT_DATASET="mock5"
DEFAULT_DATA_DIR="/gpfs/work2/0/gs/danique/mockrobiota/mock5/fastq_demux/primer"
DEFAULT_FORWARD_PRIMER="GTGCCAGCMGCCGCGGTAA"
DEFAULT_REVERSE_PRIMER="GGACTACHVGGGTWTCTAAT"
DEFAULT_REVCOMP="ATTAGAWACCCBDGTAGTCC"
DEFAULT_PRIMER_ATTACHED=1
DEFAULT_TRUNC_F=190
DEFAULT_TRUNC_R=135
DEFAULT_REMOVEOLD=0

while getopts ":d:i:f:r:p:a:t:u:y:h:" opt; do
  case $opt in
    d) DATASET="$OPTARG";;
    i) DATA_DIR="$OPTARG";;
    f) FORWARD_PRIMER="$OPTARG";;
    r) REVERSE_PRIMER="$OPTARG";;
    p) REVCOMP="$OPTARG";;
    a) PRIMER_ATTACHED="$OPTARG";;
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
      echo "  -p <reverse_complement> Specify the reverse complement of the reverse primer"
      echo "  -a <primer_attached>    Specify whether primers are still attached."
      echo "  -t <truncate_forward>   Specify at which position to truncate forward reads" 
      echo "  -u <truncate_reverse>   Specify at which position to truncate reverse reads" 
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
REVCOMP=${REVCOMP:-$DEFAULT_REVCOMP}
PRIMER_ATTACHED=${PRIMER_ATTACHED:-$DEFAULT_PRIMER_ATTACHED}
TRUNC_F=${TRUNC_F:-$DEFAULT_TRUNC_F}
TRUNC_R=${TRUNC_R:-$DEFAULT_TRUNC_R}
REMOVEOLD=${REMOVEOLD:-$DEFAULT_REMOVEOLD}

if [ "$PRIMER_ATTACHED" -eq 0 ]; then
  echo "ERROR - Sequences need to have primers attached"
  exit 1
fi 

# Pipeline Directories
BASE_DIR="/gpfs/work2/0/gs/danique/pipelines/pipelines/${PIPELINE}"
METAWORKS_DIR="${BASE_DIR}/MetaWorks1.13.0"
SNAKEFILE="${BASE_DIR}/MetaWorks1.13.0/snakefile_ESV"

INTERMEDIATE_DIR="${METAWORKS_DIR}/intermediate/${DATASET}"
OUTPUT_DIR="${BASE_DIR}/output/${DATASET}"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${INTERMEDIATE_DIR}"

# Remove files from previous run if REMOVEOLD is set to 1
if [ $REMOVEOLD -eq 1 ]; then
  echo "removing intermediate files from previous run on ${DATASET}"
  rm -rf "${INTERMEDIATE_DIR}"/* # ensure dir is empty
fi

############################################################################################
#                                 Adjust configuration file                                #
############################################################################################

cd $METAWORKS_DIR
mkdir -p "${INTERMEDIATE_DIR}/DATA"

# Add symbolic link to data folder
# Metaworks keeps all parts of the sequence ID prior to ;. If this is the sequence ID it is seen as the sample identifer and parsing fails
for f in "$DATA_DIR"/*fastq*; do
    base=$(basename "$f")
    out="$INTERMEDIATE_DIR/DATA/${base%.gz}"
    echo $out

    if [[ "$f" == *.gz ]]; then
        # stream decompress → sed → recompress
        gzip -cd "$f" | sed '1~4s/^@/@V4;/' > "$out"
        gzip -c $out > "${out}.gz"
        rm $out
    else
        sed '1~4s/^@/@V4;/' "$f" > "$out"
    fi
done

# Create a copy of the yaml file to adjust
CONFIG_FILE="${INTERMEDIATE_DIR}/${DATASET}.yaml"
cp "${BASE_DIR}/MetaWorks1.13.0/config_ESV.yaml" "${CONFIG_FILE}"

NEW_DATA_DIR="intermediate/${DATASET}/DATA"
NEW_INT_DIR="intermediate/${DATASET}/SSU"

yq --arg NEW_DATA_DIR "${NEW_DATA_DIR}" '.raw = $NEW_DATA_DIR' -i -y ${CONFIG_FILE} # Data directory
yq --arg NEW_DATA_DIR "${NEW_DATA_DIR}" '.raw_sample_read_wildcards = $NEW_DATA_DIR + "/{sample}_L001_R{read}_001.fastq.gz"' -i -y ${CONFIG_FILE} # sample format
yq --arg NEW_DATA_DIR "${NEW_DATA_DIR}" '.raw_sample_forward_wildcard = $NEW_DATA_DIR + "/{sample}_L001_R1_001.fastq.gz"' -i -y ${CONFIG_FILE} # forward format
yq --arg NEW_DATA_DIR "${NEW_DATA_DIR}" '.raw_sample_reverse_wildcard = $NEW_DATA_DIR + "/{sample}_L001_R2_001.fastq.gz"' -i -y ${CONFIG_FILE} # reverse format
yq --arg NEW_INT_DIR "${NEW_INT_DIR}" '.dir = $NEW_INT_DIR' -i -y ${CONFIG_FILE} # name of output dir
yq --arg THREADS "$THREADS" '.VSEARCH_TABLE.t = ($THREADS | tonumber)' -i -y ${CONFIG_FILE} # number of threads to use 
yq '.marker = "16S"' -i -y ${CONFIG_FILE} # name of output dir
yq '.RDP.custom = "no"' -i -y ${CONFIG_FILE} # name of output dir
yq '.pseudogene_filtering = "no"' -i -y ${CONFIG_FILE} # name of output dir

# Primer file
PRIMER_FILE="intermediate/${DATASET}/primers.fasta"
echo -e ">16sV4\n^${FORWARD_PRIMER}...${REVCOMP}$" > $PRIMER_FILE
yq --arg PRIMER_FILE "${PRIMER_FILE}" '.CUTADAPT.fasta = $PRIMER_FILE' -i -y ${CONFIG_FILE} # name of primer file

############################################################################################
#                                      Run Pipeline                                        #
############################################################################################

conda activate MetaWorks_v1.13.0
snakemake --jobs ${THREADS} --snakefile ${SNAKEFILE} --configfile ${CONFIG_FILE} --unlock
snakemake --jobs ${THREADS} --snakefile ${SNAKEFILE} --configfile ${CONFIG_FILE} --delete-all-output --rerun-incomplete
snakemake --jobs ${THREADS} --snakefile ${SNAKEFILE} --configfile ${CONFIG_FILE}
conda deactivate

############################################################################################
#                                      Export output                                       #
############################################################################################

# rm $INTERMEDIATE_DIR/DATA/*fastq*

INTERMEDIATE_DIR="${BASE_DIR}/intermediate/${DATASET}"
mkdir -p ${INTERMEDIATE_DIR}
ln -sf ${NEW_INT_DIR}  ${INTERMEDIATE_DIR}

# dna-sequences
awk -F, 'NR>1 {print ">" $1 "\n" $4}' ${NEW_INT_DIR}/results.csv > ${NEW_INT_DIR}/dna-sequences.fasta
awk 'NR%2==1 {id=$0} NR%2==0 {if (!seen[id]++) print id "\n" $0}' ${NEW_INT_DIR}/dna-sequences.fasta > ${OUTPUT_DIR}/dna-sequences.fasta

# feature table
awk -F',' 'NR>1 {
    gsub("_16s[^,]*", "", $2);  # Remove the substring from SampleName
    counts[$1][$2]=$3; 
    samples[$2]; 
    otus[$1]
} 
END {
    printf "OTU";
    for (s in samples) printf "\t%s", s;
    print "";
    for (o in otus) {
        printf "%s", o;
        for (s in samples) printf "\t%s", counts[o][s] ? counts[o][s] : "0";
        print "";
    }
}' ${NEW_INT_DIR}/results.csv > ${OUTPUT_DIR}/feature-table.tsv

# Taxonomy table
awk -F',' 'NR>1 {
    taxonomy = $9";"$12";"$15";"$18";"$21";"$24  # Concatenate taxonomic levels
    gsub(/;$/, "", taxonomy)  # Remove trailing semicolon if present
    taxa[$1] = taxonomy
} 
END {
    print "OTU\tTaxonomy"  # Print header
    for (otu in taxa) print otu "\t" taxa[otu]
}' ${NEW_INT_DIR}/results.csv > ${OUTPUT_DIR}/taxonomy.tsv

# Remove intermediate copies of the data
rm -R "${NEW_INT_DIR}/trimmed"
rm -R "${NEW_INT_DIR}/paired"

############################################################################################
#                                      Reformat output                                     #
############################################################################################

# dna-sequences
cp "${OUTPUT_DIR}/dna-sequences.fasta" "${OUTPUT_DIR}/$(basename dna-sequences.fasta .fasta)_post.fasta"

# feature table
sed 's/^OTU/ID/' "${OUTPUT_DIR}/feature-table.tsv" | \
awk 'BEGIN {FS=OFS="\t"} 
     NR==1 {for (i=1; i<=NF; i++) {split($i, a, "_"); $i=a[1]}} 
     {print}' > "${OUTPUT_DIR}/$(basename feature-table.tsv .tsv)_post.tsv" 

# taxonomy
sed 's/ /_/g' "${OUTPUT_DIR}/taxonomy.tsv" |
sed 's/;/\t/g' | \
sed '1d' | \
awk 'BEGIN {OFS="\t"; print "ID", "kingdom", "phylum", "class", "order", "family", "genus", "species"} 
     {print $1, "NA", $2, $3, $4, $5, $6, $7}' > "${OUTPUT_DIR}/$(basename taxonomy.tsv .tsv)_post.tsv" 

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
