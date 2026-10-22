#!/bin/bash
#SBATCH -p genoa
#SBATCH -t 02:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=24
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/log/sim/generatereads_%A_%a.out
#SBATCH --job-name=generatereads

# Input arguments
sample_file=$1
fasta=$2
implanted_counts=$3

if [ -z "$SLURM_CPUS_PER_TASK" ]; then
  THREADS=96  # 128 GB in MB
else
  # THREADS=$SLURM_CPUS_PER_TASK
  THREADS=${SLURM_CPUS_PER_TASK}
  echo "Threads: $THREADS"
fi

# Directories
tmp_dir="/projects/0/gs/danique/pipelines/simulations/input/tmp"
chim_dir="/projects/0/gs/danique/pipelines/simulations/input/chimeras"
fastq_dir="/projects/0/gs/danique/pipelines/simulations/input/fastq"

mkdir -p $tmp_dir
mkdir -p $chim_dir
mkdir -p $fastq_dir

# Get sample name by array index
sample=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" $sample_file)

# Subset count file for this sample
awk -v col=$((SLURM_ARRAY_TASK_ID + 2)) 'BEGIN {FS=OFS="\t"} NR>1 {print $1, $col}' "$implanted_counts" > "${tmp_dir}/${sample}_count_grinder.tsv"
total=$(awk '{sum += $2} END {print sum}' "${tmp_dir}/${sample}_count_grinder.tsv")

# Relabund
awk -v total="$total" 'BEGIN {OFS="\t"} { $2 = ($2 / total) * 100; print }' "${tmp_dir}/${sample}_count_grinder.tsv" > "${tmp_dir}/${sample}_relabund_grinder.tsv"

# First add chimeras to the fasta file (using grinder)
grinder \
  -reference_file $fasta \
  -abundance_file "${tmp_dir}/${sample}_relabund_grinder.tsv" \
  -total_reads $total \
  -read_dist 500 uniform 0 \
  -chimera_perc 20 \
  -chimera_dist 90 10 \
  -chimera_kmer 0 \
  -fastq_output 0 \
  -base_name "${sample}_chim" \
  -output_dir $chim_dir \
  -unidirectional 1 # no reverse complements 

# output file = ${sample}_chim-reads.fa
fasta_chim="${chim_dir}/${sample}_chim-simple.fa"
awk '/^>/{print ">seq" ++i; next} {print}' "${chim_dir}/${sample}_chim-reads.fa" > ${fasta_chim}

total_chim=$(grep -c '^>' $fasta_chim)

if [ ! $total_chim -eq $total ]; then
        echo "Number of sequences in chimera fasta ($total_chim) does not match number of sequences in count table ($total)."
        exit 1
fi

# Dereplicate and generate a count file
vsearch \
    --derep_fulllength "${chim_dir}/${sample}_chim-simple.fa" \
    --minuniquesize 1 \
    --sizein \
    --sizeout \
    --fasta_width 0 \
    --output "${chim_dir}/${sample}_chim-derep.fa"

# Create a count table
grep "^>" "${chim_dir}/${sample}_chim-derep.fa" \
  | sed -E 's/^>([^;]+).*size=([0-9]+).*/\1\t\2/' \
  > ${tmp_dir}/${sample}_isscount.txt

awk '{print $1 "\t" $2 * 2}' ${tmp_dir}/${sample}_isscount.txt > ${tmp_dir}/${sample}_isscount_doubled.txt

# remove size information from fasta header
sed -i 's/;size=[0-9]\+//' "${chim_dir}/${sample}_chim-derep.fa"

fasta_chim="${chim_dir}/${sample}_chim-derep.fa"

# Next generate reads using ISS with either a MiSeq or NovaSeq error model
# Activate environment
CONDA_BASE=$(conda info --base)
source $CONDA_BASE/etc/profile.d/conda.sh
conda activate iss

# Generate reads for both profiles
for profile in miseq novaseq; do

    if [ "$profile" == "miseq" ]; then
        errmod="/projects/0/gs/danique/pipelines/simulations/input/iss_mod/miseq250.npz"
    else
        errmod="novaseq"
    fi

    echo "Generating reads for $profile with error model $errmod for sample $sample"

    iss generate \
        --genomes "$fasta_chim" \
        --readcount_file ${tmp_dir}/${sample}_isscount_doubled.txt \
        --sequence_type amplicon \
        --model "$errmod" \
        --output "${fastq_dir}/${profile}/${sample}" \
        --compress \
        --cpus $THREADS

    if [ -f "${fastq_dir}/${profile}/${sample}_R1.fastq.gz" ]; then
        mv "${fastq_dir}/${profile}/${sample}_R1.fastq.gz" "${fastq_dir}/${profile}/${sample}_L001_R1_001.fastq.gz"
        mv "${fastq_dir}/${profile}/${sample}_R2.fastq.gz" "${fastq_dir}/${profile}/${sample}_L001_R2_001.fastq.gz"
        rm ${fastq_dir}/${profile}/${sample}.iss.tmp.*.vcf
    else
        echo "Error generating reads for $profile"
    fi   

    nline_fwd=$(zcat "${fastq_dir}/${profile}/${sample}_L001_R1_001.fastq.gz" | awk 'NR % 4 == 1 {print $1}' | wc -l)
    nline_rev=$(zcat "${fastq_dir}/${profile}/${sample}_L001_R2_001.fastq.gz" | awk 'NR % 4 == 1 {print $1}' | wc -l)

    if [ ! $total_chim -eq $nline_fwd ]; then
        echo "Number of sequences in chimera fasta ($total_chim) does not match reads in forward ${profile} ($nline_fwd)."
        exit 1
    fi

    if [ ! $total_chim -eq $nline_fwd ]; then
        echo "Number of sequences in chimera fasta ($total_chim) does not match reads in reverse ${profile} ($nline_rev)."
        exit 1
    fi

done

rm "${chim_dir}/${sample}"*
rm "${tmp_dir}/${sample}"*

conda deactivate