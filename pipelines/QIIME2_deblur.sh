#!/bin/bash
#SBATCH -p rome
#SBATCH -t 2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/pipelines/QIIME2Deblur/log/QIIME2deblur_output_%j.out

############################################################################################
#                              Pipeline configuration                                      #
############################################################################################

PIPELINE="QIIME2Deblur"
SILVA="https://data.qiime2.org/2023.7/common/silva-138-99-515-806-nb-classifier.qza"
SILVA_NAME=$(basename $SILVA)

############################################################################################
#                                   Run configuration                                      #
############################################################################################

if [ -z "$SLURM_CPUS_PER_TASK" ]; then
  THREADS=96  # 128 GB in MB
else
  # THREADS=$SLURM_CPUS_PER_TASK
  THREADS=$(($SLURM_NTASKS * $SLURM_CPUS_PER_TASK))
fi

# load QIIME2 software (pre-installed in snellius)
module load 2023
module load QIIME2/2023.7

# Set variables (which data to use, where to find it, and which primers are attached)
# Set default values
DEFAULT_DATASET="mock5"
DEFAULT_DATA_DIR="/gpfs/work2/0/gs/danique/mockrobiota/mock5/fastq_demux/"
DEFAULT_FORWARD_PRIMER="GTGCCAGCMGCCGCGGTAA"
DEFAULT_REVERSE_PRIMER="GGACTACHVGGGTWTCTAAT"
DEFAULT_PRIMER_ATTACHED=0
DEFAULT_REMOVEOLD=0

while getopts ":d:i:f:r:a:t:u:y:h:" opt; do
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
      echo "  -d <dataset>          Specify dataset. This will be the directory name in which intermediate (intermediate/<DATASET>) and output (output/<DATASET>) files are stored"
      echo "  -i <data_dir>         Specify data directory. This directory should contain all fastq.gz files to be processed"
      echo "  -f <forward_primer>   Specify forward primer"
      echo "  -r <reverse_primer>   Specify reverse primer"
      echo "  -a <primer_attached>  Specify whether primers are still attached."
      echo "  -y <remove_old>        Specify whether to remove old files (0 = no, 1 = yes)" 
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
REMOVEOLD=${REMOVEOLD:-$DEFAULT_REMOVEOLD}

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
  rm -rf "${INTERMEDIATE_DIR}"/* # ensure dir is empty
fi

# Download silva DB
if [ ! -f "${REF_DIR}/silva-138-99-515-806-nb-classifier.qza" ]; then
  wget -P ${REF_DIR} ${SILVA}
fi

############################################################################################
#                                      make manifest file                                  #
############################################################################################

echo "Creating manifest file..."

mkdir -p "${INTERMEDIATE_DIR}/data"
# Add symbolic link to data folder (ensuring the datafolder only contains fastq.gz files)
find "$DATA_DIR" -maxdepth 1 -name "*.fastq.gz" -exec ln -sf {} "${INTERMEDIATE_DIR}/data" \;
NEW_DATA_DIR="${INTERMEDIATE_DIR}/data"

# Create manifest file
echo -e "sample-id\tforward-absolute-filepath\treverse-absolute-filepath" > "${INTERMEDIATE_DIR}/manifest.tsv"
for f in "$NEW_DATA_DIR"/*_R1*.fastq.gz; do
  sample=$(basename "$f" | cut -d'_' -f1)
  forward="$f"
  reverse="${f/_R1/_R2}"
  echo -e "${sample}\t${forward}\t${reverse}" >> "${INTERMEDIATE_DIR}/manifest.tsv"
done

head "${INTERMEDIATE_DIR}/manifest.tsv"

############################################################################################
#                                      Run Pipeline                                        #
############################################################################################
echo "Importing fastq files..."

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "${INTERMEDIATE_DIR}/manifest.tsv" \
  --output-path "${INTERMEDIATE_DIR}/demux-paired-end.qza" \
  --input-format PairedEndFastqManifestPhred33V2

# If the primers are attached, trim them off
if [ $PRIMER_ATTACHED -eq 1 ]; then
  echo "Trimming of the primers..."

  qiime cutadapt trim-paired \
    --i-demultiplexed-sequences "${INTERMEDIATE_DIR}/demux-paired-end.qza" \
    --p-front-f "^${FORWARD_PRIMER}" \
    --p-front-r "^${REVERSE_PRIMER}" \
    --p-error-rate 0.1 \
    --o-trimmed-sequences "${INTERMEDIATE_DIR}/demux-paired-end-trimmed.qza" \
    --p-cores ${THREADS}
    
  INPUT="${INTERMEDIATE_DIR}/demux-paired-end-trimmed.qza"

  rm "${INTERMEDIATE_DIR}/demux-paired-end.qza" # remove the original demux file to save space   
else  
  INPUT="${INTERMEDIATE_DIR}/demux-paired-end.qza"
fi  

qiime vsearch merge-pairs \
  --i-demultiplexed-seqs "${INPUT}" \
  --o-merged-sequences "${INTERMEDIATE_DIR}/demux-joined.qza"

# Visualize
qiime demux summarize \
  --i-data  "${INTERMEDIATE_DIR}/demux-joined.qza" \
  --o-visualization  "${INTERMEDIATE_DIR}/demux-joined.qzv"

# Quality filter
qiime quality-filter q-score \
  --i-demux "${INTERMEDIATE_DIR}/demux-joined.qza" \
  --o-filtered-sequences "${INTERMEDIATE_DIR}/demux-joined-filtered.qza" \
  --o-filter-stats "${INTERMEDIATE_DIR}/demux-joined-filter-stats.qza"  

# Deblur
echo "Denoising using DEBLUR..."

qiime deblur denoise-16S \
  --i-demultiplexed-seqs "${INTERMEDIATE_DIR}/demux-joined-filtered.qza"  \
  --p-trim-length 252 \
  --p-sample-stats \
  --o-representative-sequences "${INTERMEDIATE_DIR}/rep-seqs-deblur.qza" \
  --o-table "${INTERMEDIATE_DIR}/table-deblur.qza" \
  --o-stats "${INTERMEDIATE_DIR}/stats-deblur.qza" \
  --p-jobs-to-start ${THREADS}

# summarize and visualize dada2 performance
qiime deblur visualize-stats \
  --i-deblur-stats "${INTERMEDIATE_DIR}/stats-deblur.qza" \
  --o-visualization "${INTERMEDIATE_DIR}/stats-deblur.qzv"

qiime feature-table summarize \
  --i-table "${INTERMEDIATE_DIR}/table-deblur.qza" \
  --o-visualization "${INTERMEDIATE_DIR}/table-deblur.qzv" \

qiime feature-table tabulate-seqs \
  --i-data "${INTERMEDIATE_DIR}/rep-seqs-deblur.qza" \
  --o-visualization "${INTERMEDIATE_DIR}/rep-seqs-deblur.qzv"

# Assign taxonomy
echo "Assigning taxonomy using a Naive Bayes Classifier (trained on SILVA)..."

qiime feature-classifier classify-sklearn \
  --i-classifier "${REF_DIR}/${SILVA_NAME}" \
  --i-reads "${INTERMEDIATE_DIR}/rep-seqs-deblur.qza" \
  --o-classification "${INTERMEDIATE_DIR}/taxonomy.qza" \
  --verbose

qiime metadata tabulate \
  --m-input-file "${INTERMEDIATE_DIR}/taxonomy.qza" \
  --o-visualization "${INTERMEDIATE_DIR}/taxonomy.qzv"

# Export files 
echo "Exporting files to tsv and biom..."

# Export the feature table and representative sequences to biom and tsv
qiime tools export \
  --input-path "${INTERMEDIATE_DIR}/table-deblur.qza" \
  --output-path "${INTERMEDIATE_DIR}" # creates a file called feature-table.biom
   
biom convert \
  -i "${INTERMEDIATE_DIR}/feature-table.biom" \
  -o "${OUTPUT_DIR}/feature-table.tsv" \
  --to-tsv

qiime tools export \
  --input-path "${INTERMEDIATE_DIR}/rep-seqs-deblur.qza" \
  --output-path "${OUTPUT_DIR}" # creates a file called dna-sequences.fasta

# Export the taxonomy table to tsv
qiime tools export \
--input-path "${INTERMEDIATE_DIR}/taxonomy.qza" \
--output-path "${INTERMEDIATE_DIR}/" # creates a file called taxonomy.tsv

# Change header of taxonomy file
sed "1s/.*/#OTUID\ttaxonomy\tconfidence/" "${INTERMEDIATE_DIR}/taxonomy.tsv" 
head "${INTERMEDIATE_DIR}/taxonomy.tsv" -n 1 | od -c

# Add taxonomy to the biom file
biom add-metadata \
-i "${INTERMEDIATE_DIR}/feature-table.biom" \
-o "${OUTPUT_DIR}/feature-table-with-tax.biom" \
--observation-metadata-fp "${INTERMEDIATE_DIR}/taxonomy.tsv" \
--observation-header ASVID,taxonomy,confidence

# Convert to tsv for easy visualization
biom convert \
-i "${OUTPUT_DIR}/feature-table-with-tax.biom" \
-o "${OUTPUT_DIR}/feature-table-with-tax.txt" \
--to-tsv \
--header-key taxonomy

# Convert to hdf5 to import in R
biom convert \
-i "${OUTPUT_DIR}/feature-table-with-tax.biom" \
-o "${OUTPUT_DIR}/feature-table-with-tax.biom-hdf5" \
--to-hdf5 \
--header-key taxonomy

############################################################################################
#                                      Export output                                       #
############################################################################################

# Reformat dna sequences and more to output folder
awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' "${OUTPUT_DIR}/dna-sequences.fasta" > tmp.fasta  && mv tmp.fasta "${OUTPUT_DIR}/dna-sequences.fasta"

# Taxonomy 
cp "${INTERMEDIATE_DIR}/taxonomy.tsv" "${OUTPUT_DIR}/taxonomy.tsv" 

############################################################################################
#                                      Reformat output                                     #
############################################################################################

# dna-sequences
cp "${OUTPUT_DIR}/dna-sequences.fasta" "${OUTPUT_DIR}/$(basename dna-sequences.fasta .fasta)_post.fasta"

# feature table
sed '1d' "${OUTPUT_DIR}/feature-table.tsv" | \
sed '1s/#OTU //g' > "${OUTPUT_DIR}/$(basename feature-table.tsv .tsv)_post.tsv"

# taxonomy
awk -F'\t' '{ $NF=""; print $0 }' OFS='\t' "${OUTPUT_DIR}/taxonomy.tsv" | \
sed 's/;/\t/g;s/.\?__//g' | \
awk 'BEGIN {OFS="\t"; print "ID", "kingdom", "phylum", "class", "order", "family", "genus", "species"} {if (NR > 1) print}' > "${OUTPUT_DIR}/$(basename taxonomy.tsv .tsv)_post.tsv"

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
