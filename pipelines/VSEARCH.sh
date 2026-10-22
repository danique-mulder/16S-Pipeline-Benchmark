#!/bin/bash
#SBATCH -p rome
#SBATCH -t 2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/pipelines/VSEARCH/log/VSEARCH_output_%j.out

############################################################################################
#                              Pipeline configuration                                      #
############################################################################################

PIPELINE="VSEARCH"
GOLD="https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.gold.bacteria.zip"
SILVA="https://www.arb-silva.de/fileadmin/silva_databases/release_138_1/Exports/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz"

############################################################################################
#                                   Run configuration                                      #
############################################################################################

if [ -z "$SLURM_CPUS_PER_TASK" ]; then
  THREADS=96  # 128 GB in MB
else
  # THREADS=$SLURM_CPUS_PER_TASK
  THREADS=$(($SLURM_NTASKS * $SLURM_CPUS_PER_TASK))
fi

# for cutadapt
export PATH=$PATH://gpfs/home4/dmulder/miniforge3/bin/

# Set variables (which data to use, where to find it, and which primers are attached)
# Set default values
DEFAULT_DATASET="mock5"
DEFAULT_DATA_DIR="/gpfs/work2/0/gs/danique/mockrobiota/mock5/fastq_demux/"
DEFAULT_FORWARD_PRIMER="GTGCCAGCMGCCGCGGTAA"
DEFAULT_REVERSE_PRIMER="GGACTACHVGGGTWTCTAAT"
DEFAULT_PRIMER_ATTACHED=0
DEFAULT_REMOVEOLD=0
DEFAULT_MINOVLEN=200
DEFAULT_MAXDIFF=15

while getopts ":d:i:f:r:a:j:k:y:h:" opt; do
  case $opt in
    d) DATASET="$OPTARG";;
    i) DATA_DIR="$OPTARG";;
    f) FORWARD_PRIMER="$OPTARG";;
    r) REVERSE_PRIMER="$OPTARG";;
    a) PRIMER_ATTACHED="$OPTARG";;
    j) MINOVLEN="$OPTARG";;
    k) MAXDIFF="$OPTARG";;
    y) REMOVEOLD="$OPTARG";;
    h)
      echo "Usage: \$0 [options]"
      echo "Options:"
      echo "  -d <dataset>            Specify dataset. This will be the directory name in which intermediate (intermediate/<DATASET>) and output (output/<DATASET>) files are stored"
      echo "  -i <data_dir>           Specify data directory. This directory should contain all fastq.gz files to be processed"
      echo "  -f <forward_primer>     Specify forward primer"
      echo "  -r <reverse_primer>     Specify reverse primer"
      echo "  -a <primer_attached>    Specify whether primers are still attached."
      echo "  -j <min_overlap_length> Specify minimum overlap length for merging reads"
      echo "  -k <max_diff>           Specify maximum number of differences allowed in the overlap region for merging reads"
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
MINOVLEN=${MINOVLEN:-$DEFAULT_MINOVLEN}
MAXDIFF=${MAXDIFF:-$DEFAULT_MAXDIFF}
REMOVEOLD=${REMOVEOLD:-$DEFAULT_REMOVEOLD}

# Pipeline Directories
BASE_DIR="/gpfs/work2/0/gs/danique/pipelines/pipelines/${PIPELINE}"
INTERMEDIATE_DIR="${BASE_DIR}/intermediate/${DATASET}"
OUTPUT_DIR="${BASE_DIR}/output/${DATASET}"
SCRIPT_DIR="/gpfs/work2/0/gs/danique/pipelines/scripts"
REF_DIR="${BASE_DIR}/intermediate"
DB="${REF_DIR}/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz_sintax.fasta"

mkdir -p "${OUTPUT_DIR}"
mkdir -p "${INTERMEDIATE_DIR}"

# Remove files from previous run if REMOVEOLD is set to 1
if [ $REMOVEOLD -eq 1 ]; then
  echo "removing intermediate files from previous run on ${DATASET}"
  rm -rf "${INTERMEDIATE_DIR}"/* # ensure dir is empty
fi

# locate software executable
VSEARCH="${BASE_DIR}/vsearch-2.28.1/bin/vsearch" # location of vsearch installation
PERL=$(which perl) # for running perl script
                            
                                                                     
############################################################################################
#                               Preparing reference database                               #
############################################################################################

#Obtaining silva gold reference (for chimera removal) 
if [ ! -e "${REF_DIR}/gold.fasta" ]; then
    if [ ! -e silva.gold.bacteria.zip ]; then
        wget -P "${REF_DIR}/" ${GOLD}
    fi

    unzip -p "${REF_DIR}/silva.gold.bacteria.zip silva.gold.align" | \
        sed -e "s/[.-]//g" > "${REF_DIR}/gold.fasta"
fi

# Obtaining silva 16s reference (for taxonomic assignment)
if [ ! -e "${REF_DIR}/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz_sintax.fasta" ]; then
  if [ ! -e "${REF_DIR}/SILVA_138.1_SSURef_NR99_tax_silva.fasta.gz" ]; then
    wget -P "${REF_DIR}/" ${SILVA}
  fi
  
  sh $SCRIPT_DIR/silva2sintax.sh -d ${DB} 
fi

############################################################################################
#                                   Remove primers (if attached)                           #
############################################################################################

# Run cutadapt if primers are still attached
if [ $PRIMER_ATTACHED -eq 1 ]; then
  # Create the trimmed directory if it doesn't exist
  mkdir -p "${INTERMEDIATE_DIR}/trimmed"
  mkdir -p "${INTERMEDIATE_DIR}/untrimmed"

  # Loop through all R1 files
  for ffile in ${DATA_DIR}/*_R1*.fastq.gz; do
    rfile="${ffile/_R1/_R2}" # R2 file
    
    fout="${INTERMEDIATE_DIR}/trimmed/$(basename "${ffile}" .fastq.gz)_trimmed.fastq.gz"
    rout="${INTERMEDIATE_DIR}/trimmed/$(basename "${rfile}" .fastq.gz)_trimmed.fastq.gz"
    
    fout_untrimmed="${INTERMEDIATE_DIR}/untrimmed/$(basename "${ffile}" .fastq.gz)_untrimmed.fastq.gz"
    rout_untrimmed="${INTERMEDIATE_DIR}/untrimmed/$(basename "${rfile}" .fastq.gz)_untrimmed.fastq.gz"

    echo "R1: $(basename $ffile)        R2: $(basename $rfile)"
    cutadapt -g "^${FORWARD_PRIMER}" \
      -G "^${REVERSE_PRIMER}" \
      -o "${fout}" \
      -p "${rout}" \
      --untrimmed-output "${fout_untrimmed}" \
      --untrimmed-paired-output "${rout_untrimmed}" \
      -j "${THREADS}" \
      "${ffile}" "${rfile}"  

  done
  
  DATA_DIR="${INTERMEDIATE_DIR}/trimmed"
  
fi

############################################################################################
#                                      Run Pipeline                                        #
############################################################################################

# Enter subdirectory                                                            
cd $INTERMEDIATE_DIR

echo
echo Checking FASTQ format version for one file

# Get information about file type and quality scores 
$VSEARCH --fastq_chars $(ls -1 ${DATA_DIR}/*.fastq.gz | head -1)

# Process samples                                                               
for f in ${DATA_DIR}/*_R1*.fastq.gz; do

    r=$(sed -e "s/_R1/_R2/" <<< "$f")
    bs=$(basename "$f")
    s=$(cut -d_ -f1 <<< "$bs")

    echo
    echo ====================================
    echo Processing sample $s
    echo ====================================
    
    s=${s//./_}
    
    echo "new sample name: ${s}"


    $VSEARCH --fastq_mergepairs $f \
        --threads $THREADS \
        --reverse $r \
        --fastqout ${s}.merged.fastq \
        --fastq_eeout \
        --fastq_minovlen $MINOVLEN \
        --fastq_maxdiffs $MAXDIFF  
     
    echo
    echo Calculate quality statistics

    $VSEARCH --fastq_eestats ${s}.merged.fastq \
        --output ${s}.stats
    echo
    echo Quality filtering

    $VSEARCH --fastq_filter ${s}.merged.fastq \
        --fastq_maxee 0.5 \
        --fastq_minlen 225 \
        --fastq_maxlen 275 \
        --fastq_maxns 0 \
        --fastaout ${s}.filtered.fasta \
        --fasta_width 0

    echo
    echo Dereplicate at sample level and relabel with sample_n

    $VSEARCH --derep_fulllength ${s}.filtered.fasta \
        --strand plus \
        --output ${s}.derep.fasta \
        --sizein \
        --sizeout \
        --uc ${s}.derep.uc \
        --relabel $s. \
        --fasta_width 0

done

echo Sum of unique sequences in each sample: $(cat *.derep.fasta | grep -c "^>")

# At this point there should be one fasta file for each sample                  
# It should be quality filtered and dereplicated.                               

echo
echo ====================================
echo Processing all samples together
echo ====================================

echo
echo Concatenate all sequences in one file

rm -f all.derep.fasta all.nonchimeras.derep.fasta
cat *.derep.fasta > all.fasta

echo
echo Dereplicate across samples and remove singletons

$VSEARCH --derep_fulllength all.fasta \
    --minuniquesize 2 \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --uc all.derep.uc \
    --output all.derep.fasta

echo Unique non-singleton sequences: $(grep -c "^>" all.derep.fasta)


# ASV
# Note: this part is not part of the published VSEARCH pipeline/tutorial, but based on the user manual
echo
echo Infer ASVs using UNOISE

$VSEARCH --cluster_unoise all.derep.fasta \
         --unoise_alpha 2 \
         --minsize 8 \
         --sizein \
         --sizeout \
         --fasta_width 0 \
         --relabel ASV_ \
         --centroids all.asvs.fasta \
         --uc all.asvs.uc

$VSEARCH --sortbysize all.asvs.fasta \
         --output all.asvs_sorted.fasta
                
echo
echo De novo chimera detection

$VSEARCH --uchime3_denovo all.asvs_sorted.fasta \
         --sizein \
         --sizeout \
         --nonchimeras all.asvs_nochim.fasta \
         --chimeras chim.fasta
         
#Removing any masking (convert any lower to upper case)
awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' all.asvs_nochim.fasta > all.asvs_nochim_nomask.fasta  

echo
echo Map reads back to ASVs

# Adapted from https://www.nicholas-ollberding.com/post/denoising-16s-sequence-reads-using-unoise3-via-vsearch-and-qiime2/    
# This deviates from the UNOISE paper (Edgar 2016), which states that reads are mapped back onto the centroids using the same skew criterium that was used to create the centroids (but w/o creating new centroids) >> same procedure as other VSEARCH-based workflows
$VSEARCH --usearch_global all.fasta \
         --db all.asvs_nochim_nomask.fasta \
         --id 0.97 \
         --otutabout asvtab.txt \
         --biomout asvtab.biom \
         --uc asvtab.uc \
         --threads $THREADS      
         
echo "ASV pipeline completed successfully."      
    
# Taxonomic Assignment (not part of tutorial, based on user manual)
$VSEARCH --sintax all.asvs_nochim.fasta \
         --db "${DB}" \
         --randseed 1234 \
         --sintax_cutoff 0.9 \
         --sintax_random \
         --tabbedout sintax_taxonomy.tsv \
         --threads 1 # for replicability (icm with randseed)        
         
echo
echo Done

date

############################################################################################
#                                      Export output                                       #
############################################################################################

if [ $PRIMER_ATTACHED -eq 1 ]; then
  # remove files with trimmed sequences
  if [ -d "${INTERMEDIATE_DIR}/trimmed" ]; then
    rm -rf "${INTERMEDIATE_DIR}/trimmed"
  fi

  # remove files with untrimmed sequences
  if [ -d "${INTERMEDIATE_DIR}/untrimmed" ]; then
    rm -rf "${INTERMEDIATE_DIR}/untrimmed"
  fi
fi

# remove large intermediate files
rm -f "${INTERMEDIATE_DIR}/"*.{filtered.fasta,merged.fastq,derep.fasta,derep.uc}

# Export final files
awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' all.asvs_nochim.fasta > "${OUTPUT_DIR}/dna-sequences.fasta"
cp asvtab.txt "${OUTPUT_DIR}/feature-table.tsv"
cp sintax_taxonomy.tsv "${OUTPUT_DIR}/taxonomy.tsv" 

############################################################################################
#                                      Reformat output                                     #
############################################################################################

# dna-sequences (asv)
sed 's/^>\([^;]*\);.*/>\1/' "${OUTPUT_DIR}/dna-sequences.fasta" > "${OUTPUT_DIR}/$(basename dna-sequences.fasta .fasta)_post.fasta"

# feature table 
sed 's/^#OTU //' "${OUTPUT_DIR}/feature-table.tsv" > "${OUTPUT_DIR}/$(basename feature-table.tsv .tsv)_post.tsv"

# taxonomy 
awk -F'\t' '{split($2, a, " +"); print $1, a[1]}' "${OUTPUT_DIR}/taxonomy.tsv" | \
sed 's/;[^ ]*//g; s/[dpcofg]://g; s/([^)]*)//g; s/,/ /g' | tr -s ' ' '\t' | \
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