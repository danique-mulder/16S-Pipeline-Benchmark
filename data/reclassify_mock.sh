
#!/bin/bash
#SBATCH -p rome
#SBATCH -t 3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=/projects/0/gs/danique/pipelines/log/blast_output_%j.out

#  BLASTs mock reference sequences SILVA v138 to update the mock reference taxonomy.

##############################################################################################
#                                       Make blast databases                                 #
##############################################################################################

# mock 5 & 12 
for dataset in mock5 mock12; do 
    DB_DIR="/projects/0/gs/danique/pipelines/db/${dataset}"
    output_dir="/projects/0/gs/danique/mockrobiota/${dataset}/info"

    if [ ! -f "${DB_DIR}/${dataset}.ndb" ]; then
        QUERY="/projects/0/gs/danique/mockrobiota/${dataset}/info/expected-sequences.fasta"
        sed 's/ /_/g' $QUERY > "${output_dir}/expected-sequences_nospace.fasta"
        QUERY="${output_dir}/expected-sequences_nospace.fasta"
        mkdir -p $DB_DIR
        makeblastdb -in $QUERY -dbtype nucl -out "${DB_DIR}/${dataset}"
    fi

    # Cut V4 from the expected-sequences (515F and reverse complement of 806R)
    if [ ! -f "${DB_DIR}_v4/${dataset}_v4.ndb" ]; then
        QUERY="/projects/0/gs/danique/mockrobiota/${dataset}/info/expected-sequences_nospace.fasta"
        V4_QUERY="/projects/0/gs/danique/mockrobiota/${dataset}/info/expected-sequences_v4.fasta"
        MISSING_FW="/projects/0/gs/danique/mockrobiota/${dataset}/info/missing_v4_FW.txt"
        MISSING_RV="/projects/0/gs/danique/mockrobiota/${dataset}/info/missing_v4_RV.txt"
        cutadapt -e 2 -g GTGCCAGCMGCCGCGGTAA -o $V4_QUERY $QUERY 
        cutadapt -e 2 -g GTGCCAGCMGCCGCGGTAA -o $MISSING_FW $QUERY --discard-trimmed

        cutadapt -e 2 -a ATTAGAWACCCBDGTAGTCC -o $V4_QUERY $V4_QUERY
        cutadapt -e 2 -g GTGCCAGCMGCCGCGGTAA -o $MISSING_RV $QUERY --discard-trimmed

        mkdir -p ${DB_DIR}_v4
        makeblastdb -in $V4_QUERY -dbtype nucl -out "${DB_DIR}_v4/${dataset}_v4"
    fi
done

# Make blastdb from SILVA (full length) and silva V4 (515F-806R) sequences
REF_DIR="/projects/0/gs/reference_db/16s"
DB_DIR="/projects/0/gs/danique/pipelines/db/silva"

if [ ! -f "${DB_DIR}/silva.ndb" ]; then
    
    # If database not downloaded, download 
    if [ ! -f "${REF_DIR}/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz" ]; then
        wget -P ${REF_DIR} https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
    fi

    gunzip -k ${REF_DIR}/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
    
    # fasta to convert into blastdb
    QUERY="${REF_DIR}/SILVA_138.2_SSURef_NR99_tax_silva.fasta"

    mkdir -p $DB_DIR
    makeblastdb -in $QUERY -dbtype nucl -out "${DB_DIR}/silva"

    # remove fasta file (keep gzipped version)
    rm "${REF_DIR}/SILVA_138.2_SSURef_NR99_tax_silva.fasta"
fi

# ##############################################################################################
# #        Run expected sequences against SILVA 138 to find relabeled sequences                #
# ##############################################################################################
# Mock only, because the simulated sequences were selected in such a way that the expected sequences do not have alternative classifications
# First run against the full length sequence.
for dataset in mock5 mock12; do 
    
    dataset="${dataset}"
    output_dir="/projects/0/gs/danique/mockrobiota/${dataset}/info"
    silva="/projects/0/gs/danique/pipelines/db/silva/silva" # full length silva blastdb
    seq=${output_dir}/expected-sequences.fasta

    if [ -f ${output_dir}/silva138-expected-taxonomy.tsv ]; then
        echo "File ${output_dir}/silva138-expected-taxonomy.tsv already exists. Skipping..."
        continue
    fi

    # Remove spaces from sequence headers
    sed 's/ /_/g' $seq > ${TMPDIR}/tmp.fasta

    # Align the query sequences against the expected sequences
    blastn -query ${TMPDIR}/tmp.fasta -db $silva -max_target_seqs 10 -out ${output_dir}/silvablast_results.txt -outfmt 6 -num_threads 32

    temp_file="${output_dir}/silvablast_results_with_tax.txt"
    silva="/projects/0/gs/reference_db/16s/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz"

    zcat $silva | sed -n 's/^>\([^ ]*\) \(.*\)/\1\t\2/p' > ${output_dir}/silva_tax_map.tsv 

    awk 'BEGIN {FS=OFS="\t"} NR==FNR {tax[$1]=$0; next} $2 in tax {print $0, tax[$2]}' ${output_dir}/silva_tax_map.tsv ${output_dir}/silvablast_results.txt > ${output_dir}/silvablast_results_with_tax.txt
    
    awk -F'\t' '{
    seq_id = $1;
    tax_field = $14;  # Assuming taxonomy is in the 14th column
    split(tax_field, taxonomy, ";");
    if (length(taxonomy) >= 7) {
        genus = taxonomy[7];
        # Clean up genus name
        gsub(/^[ \t]+|[ \t]+$/, "", genus);
        if (genus != "") {
            pair = seq_id "\t" genus;
            if (!seen[pair]++) {
                print pair;
            }
        }
    }
}' "${output_dir}/silvablast_results_with_tax.txt"  | sort  > ${output_dir}/silva138-expected-taxonomy.tsv

done

# post-edit based silva138-expected-taxonomy.tsv files based on relabeling/reclassifications in databases
# Mock 5 
# Escherichia coli > Escherichia-Shigella (coli)
# Propionibacterium acnes > Cutibacterium acnes

# Mock 12
# Clostridium phytofermentaris > Lachnoclostridium pytofermentaris
# Clostridium xylanoforans > Anaerocolumna xylanovorans 
