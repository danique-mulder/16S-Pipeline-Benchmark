#!/bin/bash
#SBATCH -p rome
#SBATCH -t 10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/pipelines/Natrix/log/Natrix_output_%j.out

############################################################################################
#                              Pipeline configuration                                      #
############################################################################################

PIPELINE="Natrix"

############################################################################################
#                                   Run configuration                                      #
############################################################################################

if [ -z "$SLURM_CPUS_PER_TASK" ]; then
  THREADS=96 
else
  # THREADS=$SLURM_CPUS_PER_TASK
  THREADS=$(($SLURM_NTASKS * $SLURM_CPUS_PER_TASK))
fi

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
DEFAULT_CLUSTER="asv"
DEFAULT_REMOVEOLD=0
DEFAULT_WAIT=0
DEFAULT_THREADS=96

while getopts ":d:i:f:r:a:c:y:w:g:h:" opt; do
  case $opt in
    d) DATASET="$OPTARG";;
    i) DATA_DIR="$OPTARG";;
    f) FORWARD_PRIMER="$OPTARG";;
    r) REVERSE_PRIMER="$OPTARG";;
    a) PRIMER_ATTACHED="$OPTARG";;
    c) CLUSTER="$OPTARG";;
    y) REMOVEOLD="$OPTARG";;
    w) WAIT="$OPTARG";;
    g) THREADS="$OPTARG";;
    h)
      echo "Usage: \$0 [options]"
      echo "Options:"
      echo "  -d <dataset>            Specify dataset. This will be the directory name in which intermediate (intermediate/<DATASET>) and output (output/<DATASET>) files are stored"
      echo "  -i <data_dir>           Specify data directory. This directory should contain all fastq.gz files to be processed"
      echo "  -f <forward_primer>     Specify forward primer"
      echo "  -r <reverse_primer>     Specify reverse primer"
      echo "  -a <primer_attached>    Specify whether primers are still attached."
      echo "  -c <cluster>         Specify how to cluster the sequences (asv, otu)."
      echo "  -y <remove_old>         Remove old files (1) or keep them (0). Default is 0" 
      echo "  -w <wait>               Wait 3 minutes before running the pipeline. Default is 0"
      echo "  -g <threads>            Specify the number of threads to use. Default is 96"
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
CLUSTER=${CLUSTER:-$DEFAULT_CLUSTER}
REMOVEOLD=${REMOVEOLD:-$DEFAULT_REMOVEOLD}
WAIT=${WAIT:-$DEFAULT_WAIT}
THREADS=${THREADS:-$DEFAULT_THREADS}

# Pipeline Directories
SCRIPT_DIR="/gpfs/work2/0/gs/danique/pipelines/scripts/"
BASE_DIR="/gpfs/work2/0/gs/danique/pipelines/pipelines/${PIPELINE}"
INTERMEDIATE_DIR="${BASE_DIR}/intermediate/${DATASET}/${CLUSTER}"
OUTPUT_DIR="${BASE_DIR}/output/${DATASET}"

mkdir -p ${INTERMEDIATE_DIR}
mkdir -p ${OUTPUT_DIR}

if [ $REMOVEOLD -eq 1 ]; then
  echo "removing intermediate files from previous run on ${DATASET}"
  rm -rf "${INTERMEDIATE_DIR}"/* # ensure dir is empty
fi

############################################################################################
#                                       Link Natrix                                        #  
############################################################################################

# If WAIT is set to 1 (pause before running for avoid conflits with other occurrences of Natrix), wait 1 minute 
if [ $WAIT -eq 1 ]; then
  sleep 180
fi

# Copy required files to the temporary directory
# This is necessary to avoid overwriting files if multiple Natrix runs are ongoing

ln -s "${BASE_DIR}/Natrix/"* "$INTERMEDIATE_DIR/"
ln -s "${BASE_DIR}/Natrix/.snakemake" "$INTERMEDIATE_DIR/"

############################################################################################
#                             Link data and create primer file                             #
############################################################################################

# Files should be named: samplename_A_R[1-2].fastq.gz
NEW_DATA_DIR="${INTERMEDIATE_DIR}/${DATASET}"
mkdir -p ${NEW_DATA_DIR}

echo "Probe,poly_N,Barcode_forward,specific_forward_primer,poly_N_rev,Barcode_reverse,specific_reverse_primer" > "${INTERMEDIATE_DIR}/${DATASET}.csv"

for direction in R1 R2; do
  for pathname in ${DATA_DIR}/*_${direction}*.fastq.gz; do
    filename=$(basename "$pathname") # extract filename
    sample=$(echo "$filename" | awk -F'_' '{print $1}' | awk -F'.' '{print $NF}') # extract sample name
    sample="${sample}_A"
    
    new_filename="${sample}_${direction}.fastq.gz"
    echo "old filename: ${filename}"
    echo "new filename: ${new_filename}"
    echo "save to:      ${NEW_DATA_DIR}/${new_filename}"
    echo
    
    cp $pathname "$NEW_DATA_DIR/${new_filename}"
    
    if [ $direction == R1 ]; then
      if [ $PRIMER_ATTACHED -eq 1 ]; then
        echo "${sample},,,${FORWARD_PRIMER},,,${REVERSE_PRIMER}" >> "${INTERMEDIATE_DIR}/${DATASET}.csv"
      else
        echo "${sample},,,,,," >> "${INTERMEDIATE_DIR}/${DATASET}.csv"
      fi  
    fi  
  done
done

head "${INTERMEDIATE_DIR}/${DATASET}.csv" # check primer file

# Change sequence headers for simulation datasets in OTU setting
max_jobs=20
job_count=0

if [[ ${DATASET} == *v4_sim* ]] && [ ${CLUSTER} == "otu" ]; then
  for d in ${NEW_DATA_DIR}/*R1.fastq.gz; do
    input_R1=${d}
    input_R2=${d/_R1/_R2}
    output_prefix=${d/_A_R1.fastq.gz/_B}

    echo $input_R1
    echo $input_R2
    echo $output_prefix
 
    (
      ${SCRIPT_DIR}/helper/make_illumina_headers.py $input_R1 $input_R2 $output_prefix

      mv "${output_prefix}_R1.fastq.gz" ${input_R1}
      mv "${output_prefix}_R2.fastq.gz" ${input_R2}
    ) &

    ((job_count++))
    if (( job_count >= max_jobs )); then
      wait  # wait for all background jobs before starting more
      job_count=0
    fi

  done

  wait

fi


############################################################################################
#                                 Adjust configuration file                                #
############################################################################################

CONFIG_FILE="${INTERMEDIATE_DIR}/${DATASET}.yaml"
cp "${BASE_DIR}/ASV.yaml" $CONFIG_FILE

# Data specification
yq --arg DATASET "$DATASET" '.general.filename = $DATASET' -i -y ${CONFIG_FILE}
yq --arg DATASET "$DATASET" '.general.primertable = $DATASET + ".csv"' -i -y ${CONFIG_FILE}
yq --arg DATASET "$DATASET" '.general.units = $DATASET + "units.tsv"' -i -y ${CONFIG_FILE}
yq --arg THREADS "$THREADS" '.general.cores = ($THREADS | tonumber)' -i -y ${CONFIG_FILE}
yq '.general.multiqc = false' -i -y $CONFIG_FILE

if [[ "$DATASET" != *_o* ]]; then
  echo "disabling cutoff (merge.cutoff = -1)"
  yq '.merge.cutoff = -1' -i -y $CONFIG_FILE # had to change (3 into -1), because otherwise all samples should have a count of 3 (which results in an empty table)
else
  echo "enabling cutoff (merge.cutoff = 3)"
  yq '.merge.cutoff = 3' -i -y $CONFIG_FILE 
fi

# Primers
if [ "$PRIMER_ATTACHED" -eq 1 ]; then
  yq '.qc.all_primer = false' -i -y ${CONFIG_FILE} # 'all_primer = true' indicated the sequence is free of adapter sequences
fi

if [ "$CLUSTER" == "otu" ]; then
  yq '.general.seq_rep = "OTU"' -i -y $CONFIG_FILE
fi

# Create unit.tsv file
cd $INTERMEDIATE_DIR
conda activate natrix

# Create units.tsv file
python create_dataframe.py "$CONFIG_FILE"
mv "${INTERMEDIATE_DIR}/units.tsv" "${INTERMEDIATE_DIR}/${DATASET}units.tsv"
head "${INTERMEDIATE_DIR}/${DATASET}units.tsv"

############################################################################################
#                                      Run Pipeline                                        #
############################################################################################

echo "running Natrix"

snakemake --use-conda --cores $THREADS --configfile $CONFIG_FILE --unlock
snakemake --use-conda --cores $THREADS --configfile $CONFIG_FILE --rerun-incomplete --jobs 30 # --cluster "sbatch -p rome -t 2-00:00:00 -c 128" 

conda deactivate

############################################################################################
#                                      Export output                                       #
############################################################################################

# Move output files to output folder
if [ ${CLUSTER} == asv ]; then
  awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' "${INTERMEDIATE_DIR}/results/finalData/filtered.fasta" > "${OUTPUT_DIR}/dna-sequences.fasta"
  
  tr ',' '\t' < "${INTERMEDIATE_DIR}/results/finalData/filtered_table.csv" > "${OUTPUT_DIR}/feature-table.tsv"
  tr ',' '\t' < "${INTERMEDIATE_DIR}/results/finalData/filtered_blast_table.csv" > "${OUTPUT_DIR}/taxonomy.tsv"
else
  awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' "${INTERMEDIATE_DIR}/results/finalData/representatives.fasta" > "${OUTPUT_DIR}/dna-sequences_${CLUSTER}.fasta"
  tr ',' '\t' < "${INTERMEDIATE_DIR}/results/finalData/filtered_blast_table.csv" > "${OUTPUT_DIR}/feature-table_${CLUSTER}.tsv"
  tr ',' '\t' < "${INTERMEDIATE_DIR}/results/finalData/filtered_blast_table.csv" > "${OUTPUT_DIR}/taxonomy_${CLUSTER}.tsv"


fi

############################################################################################
#                                      Reformat output                                     #
############################################################################################

if [ ${CLUSTER} == asv ]; then

  # dna-sequences (asv)
  sed 's/^>/>N/;s/^>\([^;]*\);.*/>\1/' "${OUTPUT_DIR}/dna-sequences.fasta" > "${OUTPUT_DIR}/$(basename dna-sequences.fasta .fasta)_post.fasta"
  
  # feature table (asv)s/;/\t/g
  sed '1s/^seqid/ID/;s/^>/N/g' "${OUTPUT_DIR}/feature-table.tsv" | \
  awk 'BEGIN {FS=OFS="\t"} {split($1, a, ";"); $1=a[1]; print}' | \
  cut -f1,3- | \
  awk 'BEGIN {FS=OFS="\t"} 
     NR==1 {for (i=1; i<=NF; i++) {split($i, a, "_"); $i=a[1]}} 
     {print}' > "${OUTPUT_DIR}/$(basename feature-table.tsv .tsv)_post.tsv" 

  # Taxonomy (asv)
  awk -F'\t' '
  NR==1 {
      for (i=1; i<=NF; i++) header[$i] = i;
      print "seqid", "taxonomy";  # New header
      next
  }
  {
      print $header["seqid"], $header["taxonomy"];
  }' OFS='\t' "${OUTPUT_DIR}/taxonomy.tsv" | \
  awk 'BEGIN {FS=OFS="\t"} {split($1, a, "_"); $1=a[1]; print}' | \
  sed 's/;/\t/g' | \
  awk 'BEGIN {OFS="\t"; print "ID", "kingdom", "phylum", "class", "order", "family", "genus", "species"} {if (NR > 1) print}' > "${OUTPUT_DIR}/$(basename taxonomy.tsv .tsv)_post.tsv" 

else

  # dna-sequences (otu)
  sed 's/^>/>N/;s/^>\([^;]*\);.*/>\1/' "${OUTPUT_DIR}/dna-sequences_otu.fasta" > "${OUTPUT_DIR}/$(basename dna-sequences_otu.fasta .fasta)_post.fasta"
  
  # feature table (otu)s/;/\t/g
  sed '1s/^seqid/ID/;s/^>/N/g' "${OUTPUT_DIR}/feature-table_otu.tsv" | \
  awk 'BEGIN {FS=OFS="\t"} {split($1, a, "_"); $1=a[1]; print}' | \
  cut -f1,14- | \
  awk 'BEGIN {FS=OFS="\t"} 
     NR==1 {for (i=1; i<=NF; i++) {split($i, a, "_"); $i=a[1]}} 
     {print}' > "${OUTPUT_DIR}/$(basename feature-table_otu.tsv .tsv)_post.tsv" 

  # Taxonomy (otu)
  awk -F'\t' '
  NR==1 {
      for (i=1; i<=NF; i++) header[$i] = i;
      print "seqid", "taxonomy";  # New header
      next
  }
  {
      print $header["seqid"], $header["taxonomy"];
  }' OFS='\t' "${OUTPUT_DIR}/taxonomy_otu.tsv" | \
  awk 'BEGIN {FS=OFS="\t"} {split($1, a, "_"); $1=a[1]; print}' | \
  sed 's/;/\t/g' | \
  awk 'BEGIN {OFS="\t"; print "ID", "kingdom", "phylum", "class", "order", "family", "genus", "species"} {if (NR > 1) print}' > "${OUTPUT_DIR}/$(basename taxonomy_otu.tsv .tsv)_post.tsv" 

fi

############################################################################################
#                                      Descriptive stats                                   #
############################################################################################

if [ ${CLUSTER} == asv ]; then 
  ## ASV
  # Write basic descriptive stats to txt file
  numasv=$(awk '/^>/' "${OUTPUT_DIR}/dna-sequences_post.fasta" | wc -l) 
  numgenus=$(awk -F'\t' '{print $7 }' "${OUTPUT_DIR}/taxonomy_post.tsv" | sort | uniq -c | wc -l)

  echo "asv: ${numasv}" > "${OUTPUT_DIR}/stat.txt"
  echo "genera: ${numgenus}" >> "${OUTPUT_DIR}/stat.txt"

  awk -F'\t' '{print $7 }' "${OUTPUT_DIR}/taxonomy_post.tsv" | sort | uniq -c | while read line; do
    echo "$line" >> "${OUTPUT_DIR}/stat.txt"
  done

else

  ## OTU
  # Write basic descriptive stats to txt file
  numasv=$(awk '/^>/' "${OUTPUT_DIR}/dna-sequences_otu_post.fasta" | wc -l) 
  numgenus=$(awk -F'\t' '{print $7 }' "${OUTPUT_DIR}/taxonomy_otu_post.tsv" | sort | uniq -c | wc -l)

  echo "asv: ${numasv}" > "${OUTPUT_DIR}/stat_otu.txt"
  echo "genera: ${numgenus}" >> "${OUTPUT_DIR}/stat_otu.txt"

  awk -F'\t' '{print $7 }' "${OUTPUT_DIR}/taxonomy_otu_post.tsv" | sort | uniq -c | while read line; do
    echo "$line" >> "${OUTPUT_DIR}/stat_otu.txt"
  done
  
fi


rm -R $NEW_DATA_DIR
