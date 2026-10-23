#!/bin/bash
#SBATCH -p rome
#SBATCH -t 10:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/log/taxonomy/learnTaxa_%j.out
#SBATCH --mem-per-cpu=500M

# Prepares the V4 reference sequences and taxonomies (SILVA, Greengenes2) and submits learnTaxa.R.

module load 2024
module load R/4.4.2-gfbf-2024a

# Run learn taxa
learnTaxa="/gpfs/work2/0/gs/danique/pipelines/scripts/classify/learnTaxa.R"
Rscript $learnTaxa

# # GreenGenes2 (fun once for download only)
# # Download the sequence and classifier files
# wget -P "/projects/0/gs/danique/pipelines/ref" https://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.v4.fna.qza
# wget -P "/projects/0/gs/danique/pipelines/ref" https://ftp.microbio.me/greengenes_release/2022.10/2022.10.backbone.tax.qza

# cd /projects/0/gs/danique/pipelines/ref

# # Extract FASTA
# unzip -p 2022.10.backbone.v4.fna.qza '*/data/dna-sequences.fasta' > greengenes2_sequences.fasta

# # Extract taxonomy to TSV
# unzip 2022.10.backbone.tax.qza -d nb_temp
# cp nb_temp/*/data/taxonomy.tsv greengenes2_taxonomy.tsv
# rm -rf nb_temp

# # Filter taxonomy to contain only subset in 
# # Extract sequence IDs
# seqkit seq -n greengenes2_sequences.fasta > gg2_ids.txt

# # Filter taxonomy
# head -1 greengenes2_taxonomy.tsv > greengenes2_taxonomy_v4.tsv
# awk 'NR==FNR{ids[$1]; next} $1 in ids' gg2_ids.txt greengenes2_taxonomy.tsv >> greengenes2_taxonomy_v4.tsv

# # Double check for identical seq ids 
# grep "^>" greengenes2_sequences.fasta | sed 's/^>//' | sed 's/ .*//' | sort > v4_ids_sorted.txt
# cut -f1 greengenes2_taxonomy_v4.tsv | sort > tax_ids_sorted.txt
# comm -3 v4_ids_sorted.txt tax_ids_sorted.txt

