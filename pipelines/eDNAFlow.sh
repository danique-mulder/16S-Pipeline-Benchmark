#!/bin/bash
#SBATCH -p rome
#SBATCH -t 10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/pipelines/eDNAFlow/log/eDNAFlow_output_%j.out

############################################################################################
#                              Pipeline configuration                                      #
############################################################################################

PIPELINE="eDNAFlow"
SILVA="https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"
VSEARCH="/gpfs/work2/0/gs/danique/pipelines/pipelines/VSEARCH/vsearch-2.28.1/bin/vsearch" # link to vsearch binary

############################################################################################
#                                   Run configuration                                      #
############################################################################################

if [ -z "$SLURM_CPUS_PER_TASK" ]; then
  THREADS=96
else
  # THREADS=$SLURM_CPUS_PER_TASK
  THREADS=$(($SLURM_NTASKS * $SLURM_CPUS_PER_TASK))
  echo $THREADS
fi

# Load java
module load 2023
module load Java/11

# Set variables (which data to use, where to find it, and which primers are attached)
# Set default values
DEFAULT_DATASET="mock5"
DEFAULT_DATA_DIR="/gpfs/work2/0/gs/danique/mockrobiota/mock5/fastq_demux/remux"
DEFAULT_FORWARD_PRIMER="GTGCCAGCMGCCGCGGTAA"
DEFAULT_REVERSE_PRIMER="GGACTACHVGGGTWTCTAAT"
DEFAULT_BARCODEFILE="/gpfs/work2/0/gs/danique/mockrobiota/mock5/fastq_demux/barcodes.txt"
DEFAULT_REMOVEOLD=0
DEFAULT_THREADS=96

while getopts ":d:i:f:r:b:y:g:h:" opt; do
  case $opt in
    d) DATASET="$OPTARG";;
    i) DATA_DIR="$OPTARG";;
    f) FORWARD_PRIMER="$OPTARG";;
    r) REVERSE_PRIMER="$OPTARG";;
    b) BARCODEFILE="$OPTARG";;
    y) REMOVEOLD="$OPTARG";;
    g) THREADS="$OPTARG";;
    h)
      echo "Usage: \$0 [options]"
      echo "Options:"
      echo "  -d <dataset>            Specify dataset. This will be the directory name in which intermediate (intermediate/<DATASET>) and output (output/<DATASET>) files are stored"
      echo "  -i <data_dir>           Specify data directory. This directory should contain all fastq.gz files to be processed"
      echo "  -f <forward_primer>     Specify forward primer"
      echo "  -r <reverse_primer>     Specify reverse primer"
      echo "  -b <barcode_file>       Specify directory of barcode file (path/to/barcode/file.txt)"
      echo "  -y <remove_old>         Remove old files (1) or keep them (0). Default is 0"
      echo "  -g <threads>            Specify number of threads to use. Default is 96"   
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
BARCODEFILE=${BARCODEFILE:-$DEFAULT_BARCODEFILE}
REMOVEOLD=${REMOVEOLD:-$DEFAULT_REMOVEOLD}
THREADS=${THREADS:-$DEFAULT_THREADS}

# Pipeline Directories
BASE_DIR="/gpfs/work2/0/gs/danique/pipelines/pipelines/${PIPELINE}"
INTERMEDIATE_DIR="${BASE_DIR}/intermediate/${DATASET}/"
OUTPUT_DIR="${BASE_DIR}/output/${DATASET}/"
LOG_DIR="${BASE_DIR}/log"
# REF_DIR="${BASE_DIR}/blast"
REF_DIR="${BASE_DIR}/16s"

mkdir -p ${INTERMEDIATE_DIR}
mkdir -p ${OUTPUT_DIR}

# link to executables
NEXTFLOW="${BASE_DIR}/nextflow"
EDNA_DIR="${BASE_DIR}/eDNAFlow"
EDNAFLOW="${EDNA_DIR}/eDNAFlow.nf"

echo "configuration is working"

############################################################################################
#                                   Install databases                                      #
############################################################################################

# Download reference DB
if [ ! -f "${REF_DIR}/16S_ribosomal_RNA.nsq" ]; then
  mkdir -p "${REF_DIR}"
  wget -P "${REF_DIR}/" https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz
  tar -xvzf "${REF_DIR}/16S_ribosomal_RNA.tar.gz" -C "${REF_DIR}/"
fi

############################################################################################
#                                         Link data                                        #
############################################################################################

# eDNAFLow requires multiplexed fastq files and barcode files to be in the same directory as the eDNAFlow executables
EDNA_DATA="${EDNA_DIR}/${DATASET}"
eDNAbarcodes="${EDNA_DATA}/barcodes.txt"

mkdir -p "${EDNA_DATA}"

# to avoid caching of a previous run
if [ -n "${EDNA_DATA}" ] && [ "$REMOVEOLD" -eq 1 ]; then 
  rm -rf ${EDNA_DATA}/*
fi

# copy data folder and barcode file
find "$DATA_DIR" -maxdepth 1 -name "*.fastq" -exec ln -sf {} "${EDNA_DATA}/" \;
cp ${BARCODEFILE} ${eDNAbarcodes}

############################################################################################
#                                         Run pipeline                                     #
############################################################################################

# Workflow minus taxonomic assignment
cd $EDNA_DATA

if [ $(wc -l < "$BARCODEFILE") -gt 2 ]; then
  echo "run with LULU"

  # # RUn eDNAFlow with LULU
  $NEXTFLOW run $EDNAFLOW --resume --barcode $eDNAbarcodes --primer-mismatch 1 --max-cpu ${THREADS} --blast_db "${REF_DIR}/16S_ribosomal_RNA" --lulu "${EDNA_DIR}/lulu.R" --mode 'vsearch' --vsearch $VSEARCH -profile nimbus --bindDir "${EDNA_DATA}/work /scratch-shared"

  # To use the curated zOTU table for LCA taxonomic assignment, some changes need to be made:
  # 1) Remove all occurrence of " in the curated file and 2) Add #ID\t to the first line
  sed 's/"//g' "${EDNA_DATA}/08_lulu_remux_minMatch84/curated_zotuTable.tab" | awk 'NR==1{print "#ID\t"$0; next}1' > "${EDNA_DATA}/08_lulu_remux_minMatch84/curated_zotuTable_forLCA.tab"

  # Taxonomic assignment (on LULU output)
  $NEXTFLOW run $EDNAFLOW --taxonomyAssignment --zotuTable "${EDNA_DATA}/08_lulu_remux_minMatch84/curated_zotuTable_forLCA.tab" --blastFile "${EDNA_DATA}/07_blast_remux/remux_blast_Result.tab" --lca_script "${EDNA_DIR}/LCA_taxonomyAssignment_scripts/runAssign_collapsedTaxonomy.py" -profile nimbus --bindDir "${EDNA_DATA}/work /scratch-shared"

  # export output and remove intermediate files to prevent overwriting
  mv ${EDNA_DATA}/09_taxonomyAssigned_lca_taxAssigned_results_qCov100_id97_diff1 ${EDNA_DATA}/09_taxonomyAssigned_lca_taxAssigned_results_qCov100_id97_diff1_lulu

else

  echo "run without LULU"

  # Run eDNAFlow without LULU
  $NEXTFLOW run $EDNAFLOW --resume --barcode $eDNAbarcodes --primer-mismatch 1 --max-cpu ${THREADS} --blast_db "${REF_DIR}/16S_ribosomal_RNA" --mode 'vsearch' --vsearch $VSEARCH -profile nimbus --bindDir "${EDNA_DATA}/work /scratch-shared"

fi

echo "run taxonomic assignment"
# Run taxonomic assignment on the output of eDNAFlow (no LULU)
$NEXTFLOW run $EDNAFLOW --taxonomyAssignment --zotuTable "${EDNA_DATA}/07_blast_remux/zotuTable.txt" --blastFile "${EDNA_DATA}/07_blast_remux/remux_blast_Result.tab" --lca_script "${EDNA_DIR}/LCA_taxonomyAssignment_scripts/runAssign_collapsedTaxonomy.py" -profile nimbus --bindDir /scratch-shared



############################################################################################
#                                      Export output                                       #
############################################################################################

# Remove copies of the data 
rm -rf "${EDNA_DATA}"/*.fastq

awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' "$(readlink -f ${EDNA_DATA}/06_Uniques_ZOTUs_remux/remux_zotus.fasta)" > "${OUTPUT_DIR}/dna-sequences.fasta" 
cp "$(readlink -f ${EDNA_DATA}/07_blast_remux/zotuTable.txt)" "${OUTPUT_DIR}/feature-table.tsv"
cp "$(readlink -f ${EDNA_DATA}/09_taxonomyAssigned_lca_taxAssigned_results_qCov100_id97_diff1/lca_taxAssigned_results_qCov100_id97_diff1.tab)" "${OUTPUT_DIR}/taxonomy.tsv"

if [ $(wc -l < "$BARCODEFILE") -gt 2 ]; then
  cp "$(readlink -f ${EDNA_DATA}/08_lulu_remux_minMatch84/curated_zotuTable_forLCA.tab)" "${OUTPUT_DIR}/feature-table-lulu.tsv"
  cp "$(readlink -f ${EDNA_DATA}/09_taxonomyAssigned_lca_taxAssigned_results_qCov100_id97_diff1_lulu/lca_taxAssigned_results_qCov100_id97_diff1.tab)" "${OUTPUT_DIR}/taxonomy-lulu.tsv"
fi

############################################################################################
#                                      Reformat output                                     #
############################################################################################

cp "${OUTPUT_DIR}/dna-sequences.fasta" "${OUTPUT_DIR}/$(basename dna-sequences.fasta .fasta)_post.fasta"
sed '1s/^#OTU //' "${OUTPUT_DIR}/feature-table.tsv" > "${OUTPUT_DIR}/$(basename feature-table.tsv .tsv)_post.tsv"

# Taxonomy
awk -F'\t' 'NR==1 {
  for (i=1; i<=NF; i++) header[$i] = i;
  print "ID", "kingdom", "phylum", "class", "order", "family", "genus";  # New header without species
  next
}
!seen[$header["OTU"]]++ {
  print $header["OTU"], "NA", $header["phylum"], $header["class"], $header["order"], $header["family"], $header["genus"];
}' OFS='\t' "${OUTPUT_DIR}/taxonomy.tsv" > "${OUTPUT_DIR}/$(basename taxonomy.tsv .tsv)_post.tsv"


if [ $(wc -l < "$BARCODEFILE") -gt 2 ]; then
  sed '1s/^#//' "${OUTPUT_DIR}/feature-table-lulu.tsv" > "${OUTPUT_DIR}/$(basename feature-table-lulu.tsv .tsv)_post.tsv"

  # Taxonomy (lulu)
  awk -F'\t' '
  NR==1 {
      for (i=1; i<=NF; i++) header[$i] = i;
      print "ID", "kingdom", "phylum", "class", "order", "family", "genus", "species";  # New header
      next
  }
  {
      print $header["OTU"], "NA", $header["phylum"], $header["class"], $header["order"], $header["family"], $header["genus"], $header["species"];
  }' OFS='\t' "${OUTPUT_DIR}/taxonomy-lulu.tsv" > "${OUTPUT_DIR}/$(basename taxonomy-lulu.tsv .tsv)_post.tsv"

fi


# Remove large intermediate files
if [ $REMOVEOLD -eq 1 ] && [ -f "${OUTPUT_DIR}/feature-table_post.tsv" ]; then
  for folder in "${EDNA_DATA}"/{00_,01_,02_,03_,04_,05_}*; do
    if [ -d "$folder" ]; then
      find "$folder" -type l \( -name "*fastq" -o -name "*fastq.gz" \) -exec bash -c 'rm -f "$(readlink -f "$0")"' {} \;
    fi
  done
fi

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

if [ $(wc -l < "$BARCODEFILE") -gt 2 ]; then
  ## LULU
  # Write basic descriptive stats to txt file
  numasv=$(awk '/^>/' "${OUTPUT_DIR}/dna-sequences_post.fasta" | wc -l) 
  numgenus=$(awk -F'\t' '{print $7 }' "${OUTPUT_DIR}/taxonomy-lulu_post.tsv" | sort | uniq -c | wc -l)

  echo "asv: ${numasv}" > "${OUTPUT_DIR}/stat-lulu.txt"
  echo "genera: ${numgenus}" >> "${OUTPUT_DIR}/stat-lulu.txt"

  awk -F'\t' '{print $7 }' "${OUTPUT_DIR}/taxonomy-lulu_post.tsv" | sort | uniq -c | while read line; do
    echo "$line" >> "${OUTPUT_DIR}/stat-lulu.txt"
  done
fi