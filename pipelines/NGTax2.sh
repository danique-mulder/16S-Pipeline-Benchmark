#!/bin/bash
#SBATCH -p rome
#SBATCH -t 2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/pipelines/NGTax2/log/NGTax2_output_%j.out

############################################################################################
#                              Pipeline configuration                                      #
############################################################################################

PIPELINE="NGTax2"
SILVA="https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"

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
SEPARATOR=${SEPARATOR:-$DEFAULT_SEPARATOR}
REMOVEOLD=${REMOVEOLD:-$DEFAULT_REMOVEOLD}

# Pipeline Directories
BASE_DIR="/gpfs/work2/0/gs/danique/pipelines/pipelines/${PIPELINE}"
INTERMEDIATE_DIR="${BASE_DIR}/intermediate/${DATASET}"
OUTPUT_DIR="${BASE_DIR}/output/${DATASET}"
LOG_DIR="${BASE_DIR}/log"
REF_DIR="${BASE_DIR}/intermediate"

mkdir -p ${INTERMEDIATE_DIR}
mkdir -p ${OUTPUT_DIR}

# Remove files from previous run if REMOVEOLD is set to 1
if [ $REMOVEOLD -eq 1 ]; then
  echo "removing intermediate files from previous run on ${DATASET}"
  rm -rf "${INTERMEDIATE_DIR}"/* # ensure dir is empty
fi

############################################################################################
#                                   Download databases                                     #
############################################################################################

# Make reference database
if [ ! -f "${REF_DIR}/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz" ]; then
  wget -P "$REF_DIR" ${SILVA}
fi

############################################################################################
#                                 Create mapping file                                      #
############################################################################################

echo "#sampleID\tforwardBarcodeSequence\tLibraryNumber\tDirection\tLibraryName\tProjectName\tRegion\tLocation" > "${INTERMEDIATE_DIR}/mappingfile.txt"

for filepath in ${DATA_DIR}/*R1*.fastq.gz; do 
  file=$(basename "$filepath")
  id=$(echo "$file" | awk -F"$SEPARATOR" '{print $1}' | awk -F'.' '{print $NF}') # extract sample name
  
  forward=$file
  reverse="${forward//R1/R2}"
  
  echo -e "${id}\t\t1\tp\t${forward},${reverse}\t${DATASET}\tV4\t" >> "${INTERMEDIATE_DIR}/mappingfile.txt" 
done 

head -n 3 "${INTERMEDIATE_DIR}/mappingfile.txt" # check mapping file

############################################################################################
#                                         Link data                                        #
############################################################################################

NEW_DATA_DIR="${INTERMEDIATE_DIR}/data"
mkdir -p "${NEW_DATA_DIR}"
rm -rf "${NEW_DATA_DIR}/"* # ensure this directory is empty
ln -sf ${DATA_DIR}/*.fastq.gz ${NEW_DATA_DIR} # symbolic link to data

############################################################################################
#                                      Run Pipeline                                        #
############################################################################################
cd $INTERMEDIATE_DIR

if [ $PRIMER_ATTACHED -eq 1 ]; then
  java -jar "${BASE_DIR}/NGTax-2.2.17.jar" -ngtax -folder "${NEW_DATA_DIR}/" -mapFile "${INTERMEDIATE_DIR}/mappingfile.txt" -for_p "${FORWARD_PRIMER}" -rev_p "${REVERSE_PRIMER}" -refdb "${REF_DIR}/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz" -b "${INTERMEDIATE_DIR}/${DATASET}.biom" -log "${LOG_DIR}" -prefixId asv
else
  java -jar "${BASE_DIR}/NGTax-2.2.17.jar" -ngtax -folder "${NEW_DATA_DIR}/" -mapFile "${INTERMEDIATE_DIR}/mappingfile.txt" -for_p "${FORWARD_PRIMER}" -rev_p "${REVERSE_PRIMER}" -refdb "${REF_DIR}/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz" -b "${INTERMEDIATE_DIR}/${DATASET}.biom" -log "${LOG_DIR}" -primersRemoved -prefixId asv
fi 

############################################################################################
#                                      Export output                                       #
############################################################################################

java -jar "${BASE_DIR}/NGTax-2.2.17.jar" -otu2fasta -i "${INTERMEDIATE_DIR}/${DATASET}.biom" -o "${INTERMEDIATE_DIR}/"
awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' "${INTERMEDIATE_DIR}/${DATASET}_otu_f.fasta" > "${OUTPUT_DIR}/dna-sequences-f.fasta" 
awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' "${INTERMEDIATE_DIR}/${DATASET}_otu_r.fasta" > "${OUTPUT_DIR}/dna-sequences-r.fasta" 

# feature table
biom convert -i "${INTERMEDIATE_DIR}/${DATASET}.biom" -o "${OUTPUT_DIR}/feature-table.tsv" --to-tsv

# taxonomy
biom convert -i "${INTERMEDIATE_DIR}/${DATASET}.biom" -o "${INTERMEDIATE_DIR}/taxonomy.tsv" --to-tsv --header-key taxonomy
awk -F'\t' '{print $1, $NF}' OFS='\t' "${INTERMEDIATE_DIR}/taxonomy.tsv" > "${OUTPUT_DIR}/taxonomy.tsv" 

############################################################################################
#                                      Reformat output                                     #
############################################################################################
# dna-sequences (f)
awk '{if ($0 ~ /^>/) {sub(/_.*/, "", $0)} print}' "${OUTPUT_DIR}/dna-sequences-f.fasta" > "${OUTPUT_DIR}/$(basename dna-sequences-f.fasta .fasta)_post.fasta"

# Only keep one unique sequence per asv identifyer 
seqkit rmdup "${OUTPUT_DIR}/$(basename dna-sequences-f.fasta .fasta)_post.fasta" -o "${OUTPUT_DIR}/dna-sequences_post.fasta"

# dna-sequences (r)
awk '{if ($0 ~ /^>/) {sub(/_.*/, "", $0)} print}' "${OUTPUT_DIR}/dna-sequences-r.fasta" > "${OUTPUT_DIR}/$(basename dna-sequences-r.fasta .fasta)_post.fasta"

# Only keep one unique sequence per asv identifyer 
seqkit rmdup "${OUTPUT_DIR}/$(basename dna-sequences-r.fasta .fasta)_post.fasta" -o tmp.fasta 
mv tmp.fasta "${OUTPUT_DIR}/dna-sequences-r_post.fasta"


# feature table
sed '1d;s/^#OTU //' "${OUTPUT_DIR}/feature-table.tsv" | \
awk 'BEGIN {FS=OFS="\t"} 
  NR==1 {for (i=2; i<=NF; i++) {split($i, a, "_"); $i=a[2]}} 
  {print}' > "${OUTPUT_DIR}/$(basename feature-table.tsv .tsv)_post.tsv" 

# taxonomy 
sed '1d;s/;/\t/g;s/.\?__//g' "${OUTPUT_DIR}/taxonomy.tsv" | \
awk 'BEGIN {OFS="\t"; print "ID", "kingdom", "phylum", "class", "order", "family", "genus", "species"} {if (NR > 1) print}' > "${OUTPUT_DIR}/$(basename taxonomy.tsv .tsv)_post.tsv" 

# Combine f and r dna-sequences for classification using Q2 and IDTAXA
seqtk seq -r ${OUTPUT_DIR}/dna-sequences-r_post.fasta > ${OUTPUT_DIR}/dna-sequences-r_revcomp.fasta

# Add 10 Ns at the end of every sequence in file1 (to create spacer)
awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' ${OUTPUT_DIR}/dna-sequences_post.fasta | \
awk '
  /^>/ {print; next}
  {print $0 "NNNNNNNNNN"}
' > ${OUTPUT_DIR}/dna-sequences_NNN.fasta

# Concatenate matching sequences from R1(+NNN) and R2
seqkit concat ${OUTPUT_DIR}/dna-sequences_NNN.fasta ${OUTPUT_DIR}/dna-sequences-r_revcomp.fasta | awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}'  > ${OUTPUT_DIR}/dna-sequences_concat.fasta


############################################################################################
#                                      Descriptive stats                                   #
############################################################################################

# Write basic descriptive stats to txt file
numasv=$(awk '/^>/' "${OUTPUT_DIR}/dna-sequences-f_post.fasta" | wc -l) 
numgenus=$(awk -F'\t' '{print $7 }' "${OUTPUT_DIR}/taxonomy_post.tsv" | sort | uniq -c | wc -l)

echo "asv: ${numasv}" > "${OUTPUT_DIR}/stat.txt"
echo "genera: ${numgenus}" >> "${OUTPUT_DIR}/stat.txt"

awk -F'\t' '{print $7 }' "${OUTPUT_DIR}/taxonomy_post.tsv" | sort | uniq -c | while read line; do
  echo "$line" >> "${OUTPUT_DIR}/stat.txt"
done
