#!/bin/bash
#SBATCH -p rome
#SBATCH -t 8:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/log/simulate_output_%j.out

# Parent script for the simulations. It builds the reference communities from GEMSIM controls, 
# implants group effects, generates MiSeq and NovaSeq reads, and runs QC, primer addition and the expected PICRUSt2 profile.


script_dir="/gpfs/work2/0/gs/danique/pipelines/scripts/simulation"
sim_dir="/gpfs/work2/0/gs/danique/pipelines/simulations/input"

# Input filenames
input_meta=${sim_dir}/metadata.tsv
input_feature=${sim_dir}/feature-table.tsv
input_fasta=${sim_dir}/dna-sequences.fasta

# Output filenames
output_meta=${sim_dir}/filtered-metadata.tsv
output_feature=${sim_dir}/filtered-feature-table.tsv
output_fasta=${sim_dir}/filtered-sequences.fasta

mkdir -p ${sim_dir}

#################################################################################################################
#                                               Copy files                                                      #
#################################################################################################################

echo "Simulating data"

if [ ! -f $input_meta ]; then
    echo "Metadata file not found. Exiting now..."
    exit 1
fi

cp /gpfs/work2/0/gs/danique/pipelines/pipelines/DADA2/output/gemsim/feature-table_post.tsv $input_feature
cp /gpfs/work2/0/gs/danique/pipelines/pipelines/DADA2/output/gemsim/dna-sequences_post.fasta $input_fasta

# #################################################################################################################
# #                                             Filter table                                                      #
# #################################################################################################################                                      

echo "Filter the metadata and count tables..."

# Filter samples (only control samples; group eq 0)
python ${script_dir}/filter_meta.py --metadata ${input_meta} --output ${output_meta} --column group --value C

# Filter features (prevalence > 0%, abundance > 0
python ${script_dir}/filter_features.py --feature_table ${input_feature} --metadata ${output_meta} --output ${output_feature}

# Filter sequences to only contain features that are in the filtered feature table
python ${script_dir}/filter_fasta.py --fasta ${input_fasta} --feature_table ${output_feature} --output ${output_fasta}

# #################################################################################################################
# #                                             Align with V4 Silva                                               #
# ################################################################################################################# 

# echo "Create blastn reference database..."
silva_seq=https://data.qiime2.org/2024.2/common/silva-138-99-seqs-515-806.qza
silva_seq_name=$(basename $silva_seq)
silva_tax=https://data.qiime2.org/2024.2/common/silva-138-99-tax-515-806.qza
silva_tax_name=$(basename $silva_tax)

# Download silva DB, if not yet installed
if [ ! -f "${sim_dir}/db/${silva_seq_name%.qza}.fasta" ]; then
  module load 2023
  module load QIIME2/2023.7 
 
  wget -P $sim_dir $silva_seq
  wget -P $sim_dir $silva_tax

  # Sequences (curated using Rescript: https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494)  
  qiime tools export \
    --input-path ${sim_dir}/${silva_seq_name} \
    --output-path ${sim_dir}/db

  mv "${sim_dir}/db/dna-sequences.fasta" "${sim_dir}/db/${silva_seq_name%.qza}.fasta"

  # Taxonomy
  qiime tools export \
    --input-path ${sim_dir}/${silva_tax_name} \
    --output-path ${sim_dir}/db
    
  mv "${sim_dir}/db/taxonomy.tsv" "${sim_dir}/db/${silva_tax_name%.qza}.tsv"

  rm ${sim_dir}/${silva_seq_name}
  rm ${sim_dir}/${silva_tax_name}
fi

silva_fasta="${sim_dir}/db/${silva_seq_name%.qza}.fasta"
silva_tax="${sim_dir}/db/${silva_tax_name%.qza}.tsv"

# Create blast database from unfiltered silva database (V4)
makeblastdb -in $silva_fasta -dbtype nucl -out "${sim_dir}/db/silvaV4" -parse_seqids

# Filter the silva database, so that simulation will result in unambiguous resolvable sequences
# Step 1: dereplicate sequences and log identical sequences
vsearch --derep_fulllength $silva_fasta --output "${silva_fasta%.fasta}-derep.fasta" --uc ${sim_dir}/silva_derep.uc --minuniquesize 1

# Step 2: Remove all sequences with ambiguous genus-level labels (same sequence, different genus)
python3 ${script_dir}/checkdup.py \
  --uc ${sim_dir}/silva_derep.uc \
  --fasta $silva_fasta \
  --tax ${silva_tax} \
  --output "${silva_fasta%.fasta}-filtered.fasta" \
  --tax_out "${silva_tax%.tsv}-filtered.tsv" \
  --excluded "${sim_dir}/excluded_silva.tsv"

silva_fasta_filtered="${silva_fasta%.fasta}-filtered.fasta"
silva_tax_filtered="${silva_tax%.tsv}-filtered.tsv"

# Create blast database from filtered silva database (V4)
makeblastdb -in $silva_fasta_filtered -dbtype nucl -out "${sim_dir}/db/silvaV4_filtered" -parse_seqids

echo "Align sequences against blastn reference database..."

# Step 3: Use blastn to find the closest match to the silva sequences
blastn -query "${output_fasta}" -db "$sim_dir/db/silvaV4_filtered" \
  -outfmt "6 qseqid sseqid pident length evalue bitscore" \
  -max_target_seqs 10 \
  -out "${sim_dir}/blast_results.tsv" \
  -num_threads 16

# sbatch -c 16 -t 15:00 --wrap="blastn -query /projects/0/gs/danique/pipelines/simulations/input/filtered-sequences.fasta -db /projects/0/gs/danique/pipelines/simulations/input/db/silvaV4_filtered -outfmt '6 qseqid sseqid pident length evalue bitscore' -max_target_seqs 10 -out /projects/0/gs/danique/pipelines/simulations/input/blast_results.tsv -num_threads 24"

echo "Filter blastn results and create a counttable and fasta file for simulation..."

# Filter the blast results to only contain the best hit per query sequence
python ${script_dir}/filter_blastn.py --blast "${sim_dir}/blast_results.tsv" --taxonomy "${silva_tax_filtered}" --output "${sim_dir}/blast_filtered.tsv" --count "${sim_dir}/count_table.tsv"

# Subset the reference db to only contain sequences corresponding to the filtered blast results
python ${script_dir}/subset_db.py --db "${silva_fasta_filtered}" --count "${sim_dir}/count_table.tsv" --output "${sim_dir}/tosimulate.fasta"

# Double check whether there are sequences with non-resolvable taxonomy
blastn -query "${sim_dir}/tosimulate.fasta" -db "$sim_dir/db/silvaV4" \
  -outfmt "6 qseqid sseqid pident length evalue bitscore qcovs" \
  -out "${sim_dir}/blast_alttax.tsv" \
  -perc_identity 100 \
  -num_threads 16

# sbatch -c 16 -t 01:00:00 --wrap="blastn -query /projects/0/gs/danique/pipelines/simulations/input/tosimulate.fasta -db /projects/0/gs/danique/pipelines/simulations/input/db/silvaV4 -perc_identity 100 -outfmt '6 qseqid sseqid pident length evalue bitscore qcovs' -out /projects/0/gs/danique/pipelines/simulations/input/blast_alttax.tsv -num_threads 16"

# filter for 100% query coverage
awk '$7 == 100' "${sim_dir}/blast_alttax.tsv" > "${sim_dir}/blast_alttax_qcov100.tsv"

# List all genus matches
python ${script_dir}/alternativetax.py --blast "${sim_dir}/blast_alttax.tsv" --taxonomy "${silva_tax}" --output "${sim_dir}/alternativeclassifications.tsv"

python ${script_dir}/filter_uniquegenus.py \
  --count_table "${sim_dir}/count_table.tsv" \
  --alt_class "${sim_dir}/alternativeclassifications.tsv" \
  --tosimulate "${sim_dir}/tosimulate.fasta" \
  --output_count "${sim_dir}/count_table_resolvable.tsv" \
  --output_fasta "${sim_dir}/tosimulate_resolvable.fasta"


# Merge the count table with the blast results to create a table with all sequence IDs and their corresponding count per sample
python ${script_dir}/merge_counts.py \
  --count "${sim_dir}/count_table_resolvable.tsv" \
  --feature_table ${output_feature} \
  --mapping "${sim_dir}/mapping_file.tsv" \
  --taxonomy "${sim_dir}/taxonomy.tsv" \
  --output "${sim_dir}/merged_counts.tsv"

awk -F'\t' '{split($2, parts, ";"); print $1 "\t" parts[6]}' "${sim_dir}/taxonomy.tsv" | tail -n +2 | sort > "${sim_dir}/taxonomy_genus.tsv" # Get genus level only for checking 

# Check if all IDs are present in tosimulate.fasta, merged_counts.tsv and taxonomy.tsv files
merged_ids=$(awk 'NR>1 {print $1}' "${sim_dir}/count_table_resolvable.tsv" | sort | uniq)
fasta_ids=$(grep "^>" "${sim_dir}/tosimulate_resolvable.fasta" | sed 's/>//g' | sort)
tax_ids=$(awk 'NR>1 {print $1}' "${sim_dir}/taxonomy.tsv" | sort)

# Check for any IDs in merged_counts.tsv that are not in tosimulate.fasta
missing_in_fasta=$(comm -23 <(echo "$merged_ids") <(echo "$fasta_ids") | wc -l)

# Check for any IDs in tosimulate.fasta that are not in merged_counts.tsv
missing_in_counts=$(comm -13 <(echo "$merged_ids") <(echo "$fasta_ids") | wc -l)

# Check for any IDs in merged_counts.tsv that are not in taxonomy.tsv
missing_in_tax=$(comm -13 <(echo "$merged_ids") <(echo "$tax_ids") | wc -l)

# If both missing_in_fasta and missing_in_counts are zero, print True, otherwise print False
if [[ $missing_in_fasta -eq 0 && $missing_in_counts -eq 0 && $missing_in_tax -eq 0 ]]; then
    echo "All files are complete and contain the same IDs."
else
    echo "Mismatch found between merged_counts.tsv, tosimulate.fasta and taxonomy.tsv."
    # exit 1
fi

#######################################################################################
#                        Classifying against q2-nb and idtaxa                         #
#######################################################################################

# learnTaxa on the same reference as Q2
# learntaxa="/projects/0/gs/danique/pipelines/scripts/classify/learnTaxa.sh"
# sbatch $learntaxa 

# Script for classifying the simulated reference sequences using q2-nbc and idtaxa 
classify=/projects/0/gs/danique/pipelines/scripts/simulation/classifysims.sh
# Run for q2-nb and idtaxa

sbatch -c 24 --mem-per-cpu=500M  $classify -c "idtaxa" -d "sim_v4"
sbatch -c 1 --mem=64G  $classify -c "q2" -d "sim_v4"
# Wait for the jobs to finish

# Output in /projects/0/gs/danique/pipelines/mock/tax
awk -F'\t' '{print $1"\t"$7}' /projects/0/gs/danique/pipelines/mock/tax/idtaxa_sim_v4_v4.tsv | tail -n +2 | sort > /projects/0/gs/danique/pipelines/mock/tax/idtaxa_sim_v4_v4_genus.tsv
awk -F'\t' '{print $1"\t"$7}' /projects/0/gs/danique/pipelines/mock/tax/q2_sim_v4.tsv | tail -n +2 | sort  > /projects/0/gs/danique/pipelines/mock/tax/q2_sim_v4_genus.tsv

#######################################################################################
#                            Filter incorrect classifications                         #
#######################################################################################
# Filter out sequences that could not be classified by either Q2 or idtaxa > we want fully resolvable sequences

exp="/projects/0/gs/danique/pipelines/simulations/input/taxonomy_genus.tsv"
q2="/projects/0/gs/danique/pipelines/mock/tax/q2_sim_v4_genus.tsv"
idtaxa="/projects/0/gs/danique/pipelines/mock/tax/idtaxa_sim_v4_v4_genus.tsv"

# ids are identical
paste <(cut -f1 $exp) <(cut -f1 $idtaxa) | awk '$1 != $2 { print "Mismatch on line " NR ": "$1" != "$2 }'
paste <(cut -f1 $exp) <(cut -f1 $q2) | awk '$1 != $2 { print "Mismatch on line " NR ": "$1" != "$2 }'

# Join exp with q2 on first column, output SILVA ID, exp_genus, q2_genus
join -t $'\t' -1 1 -2 1 $exp $q2 > ${sim_dir}/tmp_genus.tsv

# Join the above with idtaxa on first column, output SILVA ID, exp_genus, q2_genus, idtaxa_genus
join -t $'\t' -1 1 -2 1 ${sim_dir}/tmp_genus.tsv $idtaxa > ${sim_dir}/genus.tsv

# Add header 
(echo -e "SILVA_ID\texp\tq2\tidtaxa"; cat ${sim_dir}/genus.tsv) > ${sim_dir}/genus_header.tsv

python ${script_dir}/mergetax.py \
  --merged_tax ${sim_dir}/genus_header.tsv \
  --count "${sim_dir}/count_table_resolvable.tsv" \
  --fasta "${sim_dir}/tosimulate_resolvable.fasta" \
  --tax "${sim_dir}/taxonomy.tsv" \
  --mapping "${sim_dir}/mapping_file.tsv" \
  --out_merge "${sim_dir}/genus_accuracy.tsv" \
  --out_tax "${sim_dir}/taxonomy_v4.tsv" \
  --out_count "${sim_dir}/counts_v4.tsv" \
  --out_fasta "${sim_dir}/sequences_v4.fasta" \
  --out_mapping "${sim_dir}/mapping_v4.tsv"

sort -k5,5 ${sim_dir}/genus_accuracy.tsv > ${sim_dir}/genus_accuracy_sorted.tsv


#######################################################################################
#                          Remap ASVs to the filtered sequences                       #
#######################################################################################

# During the filtering process we discarded many ASVs that may have a perfect match with the selected silva IDs
# Remap the remaining sequences to the tosimulate_resolvable.fasta file containing all selected silva IDs 

# Get unmapped ASVs
comm -23 \
  <(cut -f1 $output_feature | sort | uniq) \
  <(cut -f2 ${sim_dir}/counts_v4.tsv | sort | uniq) \
  > ${sim_dir}/missing_ASVs.txt

seqkit grep -f ${sim_dir}/missing_ASVs.txt "${output_fasta}" > ${sim_dir}/missing_ASVs.fasta

# 1) Create blast database from filtered silva database (V4)
makeblastdb -in "${sim_dir}/sequences_v4.fasta" -dbtype nucl -out "${sim_dir}/db/tosimulate" -parse_seqids

blastn -query ${sim_dir}/missing_ASVs.fasta -db "${sim_dir}/db/tosimulate" \
  -outfmt "6 qseqid sseqid pident length evalue bitscore" \
  -out "${sim_dir}/blast_remap.tsv" \
  -perc_identity 90 \
  -num_threads 16

# sbatch -c 16 -t 01:00:00 --wrap="blastn -query /gpfs/work2/0/gs/danique/pipelines/simulations/input/missing_ASVs.fasta -db /projects/0/gs/danique/pipelines/simulations/input/db/tosimulate -perc_identity 90 -outfmt '6 qseqid sseqid pident length evalue bitscore' -out /projects/0/gs/danique/pipelines/simulations/input/blast_remap.tsv -num_threads 16"

# Combine the two files
python ${script_dir}/remap_asvs.py \
  --blast "${sim_dir}/blast_remap.tsv" \
  --taxonomy "${sim_dir}/taxonomy_v4.tsv" \
  --count "${sim_dir}/counts_v4.tsv" \
  --output "${sim_dir}/remapped_counts.tsv"

# Create a count table
python ${script_dir}/merge_counts.py \
  --count "${sim_dir}/remapped_counts.tsv" \
  --feature_table ${output_feature} \
  --mapping "${sim_dir}/mapping_file.tsv" \
  --taxonomy "${sim_dir}/taxonomy.tsv" \
  --output "${sim_dir}/merged_counts.tsv"


#######################################################################################
#                                       Final check                                   # 
#######################################################################################

# Check if all IDs are present in tosimulate.fasta, merged_counts.tsv and taxonomy.tsv files
merged_ids=$(awk 'NR>1 {print $1}' "${sim_dir}/merged_counts.tsv" | sort)
fasta_ids=$(grep "^>" "${sim_dir}/sequences_v4.fasta" | sed 's/>//g' | sort)
tax_ids=$(awk 'NR>1 {print $1}' "${sim_dir}/taxonomy.tsv" | sort)

# Check for any IDs in merged_counts.tsv that are not in tosimulate.fasta
missing_in_fasta=$(comm -23 <(echo "$merged_ids") <(echo "$fasta_ids") | wc -l)

# Check for any IDs in tosimulate.fasta that are not in merged_counts.tsv
missing_in_counts=$(comm -13 <(echo "$merged_ids") <(echo "$fasta_ids") | wc -l)

# Check for any IDs in merged_counts.tsv that are not in taxonomy.tsv
missing_in_tax=$(comm -13 <(echo "$merged_ids") <(echo "$tax_ids") | wc -l)

# If both missing_in_fasta and missing_in_counts are zero, print True, otherwise print False
if [[ $missing_in_fasta -eq 0 && $missing_in_counts -eq 0 && $missing_in_tax -eq 0 ]]; then
    echo "All files are complete and contain the same IDs."
else
    echo "Mismatch found between merged_counts.tsv, tosimulate.fasta and taxonomy.tsv."
    # exit 1
fi

#################################################################################################################
#                                                  Implant Effect                                               #
################################################################################################################# 

# Access r libraries
R_LIBS=/projects/0/gs/danique/rpackages:$R_LIBS
export R_LIBS

# Load R module (pre-installed in Snellius)
module load 2024
module load R/4.4.2-gfbf-2024a

# # Files needed for simulation
fasta="${sim_dir}/sequences_v4.fasta"
counts="${sim_dir}/merged_counts.tsv" 
tax="${sim_dir}/taxonomy.tsv"
implanted_counts="${sim_dir}/implanted_counts.tsv"

echo "Implant an effect using simbaR..."

# Perform implantation using simbaR
Rscript --vanilla ${script_dir}/simba.R \
  --count $counts \
  --metadata $output_meta \
  --taxonomy $tax \
  --abundance "1.5,2" \
  --prevalence "0.1, 0.2" \
  --feature_type "middle_high" \
  --output_count ${implanted_counts} \
  --list "${sim_dir}/sim_overview.rds" \
  --tmp "${sim_dir}" \
  --keep FALSE

# Check if all IDs are still present in both implanted_counts and tosimulate.fasta and 
implanted_ids=$(awk 'NR>1 {print $1}' ${implanted_counts} | sort)
fasta_ids=$(grep "^>" ${fasta}| sed 's/>//g' | sort)
missing_in_fasta=$(comm -23 <(echo "$implanted_ids") <(echo "$fasta_ids") | wc -l)
missing_in_counts=$(comm -13 <(echo "$implanted_ids") <(echo "$fasta_ids") | wc -l)

# If both missing_in_fasta and missing_in_counts are zero, print True, otherwise print False
if [[ $missing_in_fasta -eq 0 && $missing_in_counts -eq 0 ]]; then
    echo "All files are complete and contain the same IDs."
else
    echo "Mismatch found between merged_counts.tsv, tosimulate.fasta and taxonomy.tsv."
    # exit 1
fi  

#################################################################################################################
#                                  Simulate chimeras + simulate sequencing                                      #
################################################################################################################# 

# Install grinder
# mkdir -p /projects/0/gs/danique/biogrinder && cd /projects/0/gs/danique/biogrinder
# wget https://sourceforge.net/projects/biogrinder/files/biogrinder/Grinder-0.5.4/Grinder-0.5.4.tar.gz
# tar -xzf Grinder-0.5.4.tar.gz
# cd Grinder-0.5.4

# perl Makefile.PL PREFIX=$HOME/perl5
# make
# make install
# # Files needed for simulation

# # Install insilicoseq
# https://insilicoseq.readthedocs.io/en/latest/iss/generate.html#amplicons
# # One time setup
# conda create -n iss python=3.8 -y # create env
# conda activate iss # activate env
# conda install -c bioconda insilicoseq # Install insilicoseq

# Create a miseq 2x250 error model (only 2x300 is available, which is too long for 525-806 primer pair), if it doesn't exist already
if [ ! -f "${sim_dir}/iss_mod/miseq250.npz" ]; then
  
  echo "Build an error model..."
  
  mkdir -p "${sim_dir}/iss_mod"
  # source conda environments
  CONDA_BASE=$(conda info --base)
  CONDA_SH="$CONDA_BASE/etc/profile.d/conda.sh"
  source $CONDA_SH
  conda activate iss
   
  module load 2023
  module load Bowtie2/2.5.1-GCC-12.3.0
  module load SAMtools/1.18-GCC-12.3.0

  mock5ref=/projects/0/gs/danique/mockrobiota/mock5/info/expected-sequences_num.fasta
  mock5_R1=/projects/0/gs/danique/mockrobiota/mock5/fastq_demux/Even1_L001_R1_001.fastq.gz
  mock5_R2=/projects/0/gs/danique/mockrobiota/mock5/fastq_demux/Even1_L001_R2_001.fastq.gz

  # Align reads against the reference sequences
  bowtie2-build "$mock5ref" "${sim_dir}/iss_mod/mock5_index"
  bowtie2 -p 24 -x "${sim_dir}/iss_mod/mock5_index" -1 ${mock5_R1} -2 ${mock5_R2} \
    | samtools view -@ 24 -bS - \
    | samtools sort -@ 24 -o "${sim_dir}/iss_mod/mock5_mapping.bam"
  samtools index -@ 24 "${sim_dir}/iss_mod/mock5_mapping.bam"

  # Build an iss model
  iss model -b "${sim_dir}/iss_mod/mock5_mapping.bam" -o "${sim_dir}/iss_mod/miseq250"

  conda deactivate
fi

########################

echo "Simulate sequencing error using insilicoseq..."

# Run simulation per sample
fasta="${sim_dir}/sequences_v4.fasta"
implanted_counts="${sim_dir}/implanted_counts.tsv"
generate_reads="${script_dir}/generatereads.sh"
log_dir="/projects/0/gs/danique/pipelines/log/sim"
fastq_dir="${sim_dir}/fastq"
tmp_dir="${sim_dir}/tmp"
mkdir -p $fastq_dir
mkdir -p $fastq_dir/miseq
mkdir -p $fastq_dir/novaseq
mkdir -p $tmp_dir
mkdir -p $log_dir

# get header and split into samplenames
header=$(head -n 1 "$implanted_counts")
IFS=$'\t' read -r -a samples <<< "$header" 

for ((i=1; i<${#samples[@]}; i++)); do
  echo "${samples[i]}"
done > ${sim_dir}/sample_list.txt

# head -n 2 ${sim_dir}/sample_list.txt > tmp.txt 
# mv tmp.txt ${sim_dir}/sample_list.txt 

# Assuming the file has a header line you want to exclude
total_lines=$(wc -l < "${sim_dir}/sample_list.txt")
num_samples=$((total_lines - 1)) # 0 indexed so minus 1 

# Submit the array job with correct size
job_id=$(sbatch --array=0-${num_samples}%16 -t 30:00 -c 24 $generate_reads ${sim_dir}/sample_list.txt $fasta $implanted_counts | awk '{print $4}')


# Check sequence length for the following sequence IDs (not simulated due to length)
grep "shorter than read length for this ErrorModel" generatereads_15391458_*.out | awk -F':' '{print $4}' | awk -F' ' '{print $1}' | sort | uniq -c > "${sim_dir}/notsimulated_miseq.txt"

QC=/projects/0/gs/danique/pipelines/scripts/qc/qc.sh
sbatch -t 01:00:00 $QC -d "miseq_sim" -i "${fastq_dir}/miseq"
sbatch -t 01:00:00 $QC -d "novaseq_sim" -i "${fastq_dir}/novaseq"

# # Add primers (needed for some of the workflows)
ADDPRIMERS=/projects/0/gs/danique/pipelines/scripts/helper/addprimers.sh
sbatch -c 1 $ADDPRIMERS "${fastq_dir}/miseq"
sbatch -c 1 $ADDPRIMERS "${fastq_dir}/novaseq"

##### Get picrust-predicted functional annotations for the reference sequences 

sbatch ${script_dir}/picrust2_exp.sh

# #################################################################################################################
# #                                         Assess alternative classifications                                    #
# ################################################################################################################# 

# silva_seq=https://data.qiime2.org/2024.2/common/silva-138-99-seqs-515-806.qza
# silva_seq_name=$(basename $silva_seq)
# silva_tax=https://data.qiime2.org/2024.2/common/silva-138-99-tax-515-806.qza
# silva_tax_name=$(basename $silva_tax)

# silva_fasta="${sim_dir}/db/${silva_seq_name%.qza}.fasta"
# silva_tax="${sim_dir}/db/${silva_tax_name%.qza}.tsv"
# fasta="${sim_dir}/tosimulate.fasta"

# # Use blast database from the unfiltered silva database
# makeblastdb -in $silva_fasta -dbtype nucl -out "${sim_dir}/db/silvaV4_unfilt" -parse_seqids

# echo "Align to simulate sequences against blast reference database..."

# blastn -query "${output_fasta}" -db "${sim_dir}/db/silvaV4_unfilt" \
#   -outfmt "6 qseqid sseqid pident length evalue bitscore" \
#   -out "${sim_dir}/blast_results_simulated.tsv" \
#   -num_threads 24

# sbatch --job-name=blast_simulated \
#        --time=2:00:00 \
#        --cpus-per-task=24 \
#        --wrap="blastn -query \"${output_fasta}\" -db \"${sim_dir}/db/silvaV4_unfilt\" \
#                 -outfmt \"6 qseqid sseqid pident length evalue bitscore\" \
#                 -out \"${sim_dir}/blast_results_simulated.tsv\" \
#                 -num_threads 24"

# python ${script_dir}/alternativetax.py --blast "${sim_dir}/blast_results_simulated.tsv" --taxonomy "${silva_tax}" --output "${sim_dir}/alternativeclassifications.tsv"
