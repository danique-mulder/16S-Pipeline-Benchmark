#!/bin/bash
#SBATCH -p rome
#SBATCH -t 3:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --output=/projects/0/gs/danique/pipelines/log/blast_output_%j.out


##############################################################################################
#                                       Make blast databases                                 #
##############################################################################################

# mock 5 & 12 reference sequences
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

# create blast db for the simulated (insilicoseq; same for miseq_v4_sim and novaseq_v4_sim)
dataset=sim_v4
DB_DIR="/projects/0/gs/danique/pipelines/db/${dataset}"

if [ ! -f "${DB_DIR}/${dataset}.ndb" ]; then
    DB="/projects/0/gs/danique/pipelines/simulations/input/sequences_v4.fasta"
    mkdir -p $DB_DIR
    makeblastdb -in $DB -dbtype nucl -out "${DB_DIR}/${dataset}"
fi

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

# Same for SILVA V4 (RESCRIPT curated)
silva_seq=https://data.qiime2.org/2024.2/common/silva-138-99-seqs-515-806.qza
silva_seq_name=$(basename $silva_seq)
silva_tax=https://data.qiime2.org/2024.2/common/silva-138-99-tax-515-806.qza
silva_tax_name=$(basename $silva_tax)

# Download silva DB for V4 sequences as qza (if not already present) and convert to fasta
if [ ! -f "${REF_DIR}/${silva_seq_name%.qza}.fasta" ]; then
  module load 2023
  module load QIIME2/2023.7 
 
  wget -P $REF_DIR $silva_seq
  wget -P $REF_DIR $silva_tax

  # Sequences (curated using Rescript: https://forum.qiime2.org/t/processing-filtering-and-evaluating-the-silva-database-and-other-reference-sequence-data-with-rescript/15494)  
  qiime tools export \
    --input-path ${REF_DIR}/${silva_seq_name} \
    --output-path ${REF_DIR}

  mv "${REF_DIR}/dna-sequences.fasta" "${REF_DIR}/${silva_seq_name%.qza}.fasta"

  # Taxonomy
  qiime tools export \
    --input-path ${REF_DIR}/${silva_tax_name} \
    --output-path ${REF_DIR}/

  mv "${REF_DIR}/taxonomy.tsv" "${REF_DIR}/${silva_tax_name%.qza}.tsv"

  rm ${REF_DIR}/${silva_seq_name}
  rm ${REF_DIR}/${silva_tax_name}
fi

# Convert to blastdb
silvav4="${REF_DIR}/${silva_seq_name%.qza}.fasta"
DB_DIR="/projects/0/gs/danique/pipelines/db/silva_v4"

if [ ! -f "${DB_DIR}/silva_v4.ndb" ]; then
    makeblastdb -in $silvav4 -dbtype nucl -out "${DB_DIR}/silva_v4"
fi

##############################################################################################
#                                         Run config                                         #
##############################################################################################

r="recall precision" # recall precision
pipelines="all" # all AMPtk VSEARCH DADA2 Tourmaline dadasnake nfcore_ampliseq Natrix DADAist2 QIIME2 MetaWorks APSCALE LotuS2 NGTax2 Natrix LotuS2_otu
datasets="mock5 mock12 miseq_v4_sim novaseq_v4_sim" # mock5 mock12 miseq_v4_sim novaseq_v4_sim mock5_o mock12_o

##############################################################################################
#                                       Assess precision                                     #
##############################################################################################
# For each pipeline and each dataset, run blast against the expected sequences first. For all bad matches (>=2 mismatches), run blast against SILVA

if [[ $r == *precision* ]]; then
    echo "Assessing precision (blast observed sequences against expected sequences)..."

    silva="/projects/0/gs/danique/pipelines/db/silva_v4/silva_v4"
    output_dir="/projects/0/gs/danique/pipelines/precision"
    mkdir -p $output_dir

    for dataset in $datasets; do
        for suffix in "" "_otu"; do

            echo "Processing dataset: ${dataset}"
            
            if [[ $dataset == "miseq_v4_sim" || $dataset == "novaseq_v4_sim" ]]; then
                db="/projects/0/gs/danique/pipelines/db/sim_v4/sim_v4"
            else
                d="${dataset/_o/}"
                db="/projects/0/gs/danique/pipelines/db/${d}/${d}"
            fi


            for pipeline in $pipelines; do

                if [[ $suffix == "_otu" ]] && [[ $pipeline != "LotuS2" ]] && [[ $pipeline != "Natrix" ]]; then
                    continue
                fi

                echo "  Processing pipeline: ${pipeline}${suffix}"

                if [[ $pipeline == "all" ]]; then
                    seq="/projects/0/gs/danique/pipelines/int_results/uni_asv_${dataset}_derep.fasta"

                    intermediate_dir="/projects/0/gs/danique/pipelines/int_results/blast" 
                    mkdir -p $intermediate_dir

                else
                    seq="/projects/0/gs/danique/pipelines/pipelines/${pipeline}/output/${dataset}/dna-sequences${suffix}_post.fasta"
                    intermediate_dir="/projects/0/gs/danique/pipelines/pipelines/${pipeline}/intermediate/${dataset}/blast" 
                fi

                mkdir -p $intermediate_dir
                
                if [ ! -f $seq ]; then
                    echo "File $seq not found. Skipping..."
                    continue
                fi

                # Process files with current suffix
                if ! [ $pipeline == "NGTax2" ]; then
                    # best hit is sufficient
                    blastn -query $seq -db $db -max_target_seqs 1 -out "${intermediate_dir}/smallblast_results${suffix}.txt" -outfmt 6 -num_threads 32 

                    echo -e "ID\ttax\tmismatch\treference" > "${intermediate_dir}/smallblast_good${suffix}.tsv"

                    # the fifth column is the number of mismatches
                    awk -v OFS='\t' '
                    NR > 1 {
                        if ($5 <= 1) {
                            print $1, $2, $5, "Reference" >> "'${intermediate_dir}/smallblast_good${suffix}.tsv'"
                        } else {
                            print $1, $2, $5, "Other" > "'${intermediate_dir}/smallblast_bad${suffix}.tsv'"
                        }
                    }' "${intermediate_dir}/smallblast_results${suffix}.txt"
                else 
                    # for NGTax2, get all hits to be able to compare forward and reverse hits
                    blastn -query $seq -db $db -out "${intermediate_dir}/smallblast_results${suffix}.txt" -outfmt 6 -num_threads 32 # -max_target_seqs 1 

                    echo "Blasting reverse sequences for NG-Tax2..."

                    # For NGTax, get all IDs from both forward and reverse sequences and combine to determine bad matches
                    # First, blast the reverse sequences against the expected sequences
                    seq2="/projects/0/gs/danique/pipelines/pipelines/${pipeline}/output/${dataset}/dna-sequences-r${suffix}_post.fasta"
                    blastn -query $seq2 -db $db -strand minus -out "${intermediate_dir}/smallblast_results_r${suffix}.txt" -outfmt 6 -num_threads 32 # Only search the reverse complement (strand = minus)
                    # head "${intermediate_dir}/smallblast_results.txt"

                    # Creates smallblast_good and _bad files for further processing 
                    python /projects/0/gs/danique/pipelines/scripts/benchmark/precision_NGTax_asv.py "${intermediate_dir}/smallblast_results${suffix}.txt" "${intermediate_dir}/smallblast_results_r${suffix}.txt" "${intermediate_dir}"   
                fi

                # Create fasta file with bad sequences
                grep '^>' $seq | sed 's/^>//' > "${intermediate_dir}/all_ids${suffix}.txt"
                cut -f1 "${intermediate_dir}/smallblast_good${suffix}.tsv" | tail -n +2 > "${intermediate_dir}/good_ids${suffix}.txt"
                grep -F -x -v -f "${intermediate_dir}/good_ids${suffix}.txt" "${intermediate_dir}/all_ids${suffix}.txt" > "${intermediate_dir}/bad_ids${suffix}.txt"

                seqtk subseq $seq "${intermediate_dir}/bad_ids${suffix}.txt" > "${intermediate_dir}/bad_sequences${suffix}.fasta"

                if [ $pipeline == "NGTax2" ]; then
                    seqtk subseq $seq2 "${intermediate_dir}/bad_ids.txt" > "${intermediate_dir}/bad_sequences_r.fasta"
                fi


                # Blast the bad matcches against SILVA (V4)
                # If there are no "bad sequences", or if the sequences are simulated (no contaminants) skip the SILVA blast
                if [[ $dataset != *"_sim" ]] && [[ $pipeline != "NGTax2" ]] && [ $(wc -l < "${intermediate_dir}/bad_sequences${suffix}.fasta") -gt 0 ]; then
                    echo "Blasting against SILVA..."
                    blastn -query "${intermediate_dir}/bad_sequences${suffix}.fasta" -db $silva -max_target_seqs 1 -best_hit_score_edge 0.1 -out "${intermediate_dir}/bigblast_results${suffix}.txt" -outfmt 6 -num_threads 32
                    awk -v OFS='\t' '{print $1, $2, $5, "Silva"}' "${intermediate_dir}/bigblast_results${suffix}.txt" > "${intermediate_dir}/bigblast_all${suffix}.tsv"
                    cat "${intermediate_dir}/smallblast_good${suffix}.tsv" "${intermediate_dir}/bigblast_all${suffix}.tsv" > "${output_dir}/precision_${pipeline}_${dataset}${suffix}.tsv"    
                else
                    echo "Skipping SILVA blast..."
                    cat "${intermediate_dir}/smallblast_good${suffix}.tsv" "${intermediate_dir}/smallblast_bad${suffix}.tsv"  > "${output_dir}/precision_${pipeline}_${dataset}${suffix}.tsv"
                fi

                # Add missing query sequences to results file
                cut -f1 "${output_dir}/precision_${pipeline}_${dataset}${suffix}.tsv" | tail -n +2 > "${intermediate_dir}/present_ids${suffix}.txt"
                grep -Fxvf "${intermediate_dir}/present_ids${suffix}.txt" "${intermediate_dir}/all_ids${suffix}.txt" > "${intermediate_dir}/missing_ids${suffix}.txt"

                while read -r id; do
                    echo -e "${id}\tNA\t9999\tNoHit" >> "${output_dir}/precision_${pipeline}_${dataset}${suffix}.tsv"
                done < "${intermediate_dir}/missing_ids${suffix}.txt"
            done
        done
    done
fi

##############################################################################################
#                                       Assess recall                                        #
##############################################################################################

if [[ $r == *recall* ]]; then

    echo "Assessing recall (blast expected sequences against observed sequences)..."

    # I check for each sequence in the expected-sequences.fasta file, whether this one is represented in the dna-sequences_post.fasta file (pipeline specific)
    # Query sequences =  expected sequences
    # DB = observed sequences per pipeline 

    # Create one table per pipeline, showing for each expected sequence, with which ASV in the output it has a perfect match 
    # Rownames > expected ID

    # The output file is a table with the following columns:
    # Column 1 > ASV-ID (if multiple ; separated) with 0 mismatches
    # Column 2 > ASV-ID (if multiple ; separated) with 1 mismatch
    # Column 3 > ASV-ID (if multiple ; separated) with 2+ mismatches

    output_dir="/projects/0/gs/danique/pipelines/recall"
    mkdir -p "$output_dir"

    for dataset in $datasets; do
        for suffix in "" "_otu"; do

            echo "Processing dataset: $dataset"

            # Specify sequence file (expected sequences)    
            if [[ $dataset == *mock* ]]; then
                d="${dataset/_o/}"
                SEQ="/projects/0/gs/danique/mockrobiota/${d}/info/expected-sequences_num.fasta" # numbered, to be able to parse hits to different seq with same id
                # number the sequences
                awk '/^>/ {print ">" ++i "_" substr($0, 2); next} {print}' "/projects/0/gs/danique/mockrobiota/${d}/info/expected-sequences_nospace.fasta" > ${SEQ}
            elif [[ $dataset == "miseq_v4_sim" || $dataset == "novaseq_v4_sim" ]]; then
                SEQ="/projects/0/gs/danique/pipelines/simulations/input/sequences_v4.fasta"
            fi

            for pipeline in $pipelines; do

                if [[ $suffix == "_otu" ]] && [[ $pipeline != "LotuS2" ]] && [[ $pipeline != "Natrix" ]]; then
                    continue
                fi
                
                echo "  Processing pipeline: $pipeline${suffix}"

                # Set input/output files and directories
                DB_DIR="/projects/0/gs/danique/pipelines/db/${dataset}_${pipeline}"

                if [ $pipeline == "all" ]; then
                    DB="/projects/0/gs/danique/pipelines/int_results/uni_asv_${dataset}_derep.fasta" # dereplicated concatenated file
                    intermediate_dir="/projects/0/gs/danique/pipelines/int_results/blast"
                else 
                    DB="/projects/0/gs/danique/pipelines/pipelines/${pipeline}/output/${dataset}/dna-sequences${suffix}_post.fasta"
                    intermediate_dir="/projects/0/gs/danique/pipelines/pipelines/${pipeline}/output/${dataset}/blast"
                fi
                
                recall_file="${output_dir}/recall_${pipeline}_${dataset}${suffix}.tsv"
                
                mkdir -p $intermediate_dir
                
                # Create recall file for dataset (with the ids of all expected sequences)
                expected_ids="${intermediate_dir}/expected_ids${suffix}.txt"
                seqkit seq -n ${SEQ} > ${expected_ids}
                header="sequence-id"

                if [ ! -f $DB ]; then
                    echo "File $DB not found. Skipping..."
                    continue
                fi

                # Make tmp blast directory of the dna-sequences_post.fasta
                mkdir -p $DB_DIR
                makeblastdb -in $DB -dbtype nucl -out "${DB_DIR}/${dataset}_${pipeline}${suffix}"

                # Run blastn
                blastn \
                    -query ${SEQ} \
                    -db "${DB_DIR}/${dataset}_${pipeline}${suffix}" \
                    -outfmt "6 qseqid sseqid pident mismatch length" \
                    -out ${intermediate_dir}/recall_blast_all${suffix}.txt \
                    -num_threads 32

                if [ $pipeline == "NGTax2" ]; then
                    # Also blast the reverse sequences 
                    DB2="/projects/0/gs/danique/pipelines/pipelines/${pipeline}/output/${dataset}/dna-sequences-r${suffix}_post.fasta"
                    
                    makeblastdb -in $DB2 -dbtype nucl -out "${DB_DIR}/${dataset}_${pipeline}_r${suffix}"
                    # Run blastn
                    blastn \
                        -query ${SEQ} \
                        -db "${DB_DIR}/${dataset}_${pipeline}_r${suffix}" \
                        -outfmt "6 qseqid sseqid pident mismatch length" \
                        -out ${intermediate_dir}/recall_blast_all_r${suffix}.txt \
                        -num_threads 32
                fi


                if [ "$pipeline" != "NGTax2" ]; then
                    # Parse blast results to group ASV IDs by mismatch count category per expected sequence (0, 1, 2+ mismatches)
                    awk '
                        BEGIN { FS = "\t"; }
                        # Read expected IDs in order
                        FNR==NR { expected[$1] = 1; order[++n] = $1; next; }
                        {S
                            q = $1; s = $2; m = $4;
                            if (!(q in min_mismatch) || m < min_mismatch[q]) {
                                min_mismatch[q] = m;
                                min_asvs[q] = s;
                            } else if (m == min_mismatch[q]) {
                                min_asvs[q] = min_asvs[q] ";" s;
                            }
                            if (m == 0) {
                                zero[q] = (q in zero) ? zero[q] ";" s : s;
                            } else if (m == 1) {
                                one[q] = (q in one) ? one[q] ";" s : s;
                            }
                        }
                        END {
                            print "sequence-ID\tASV_0_mismatch\tASV_1_mismatch\tASV_min_mismatch";
                            for (i = 1; i <= n; i++) {
                                q = order[i];
                                if (q in zero) {
                                    print q "\t" zero[q] "\t\t";
                                } else if (q in one) {
                                    print q "\tNA\t" one[q] "\t";
                                } else {
                                    print q "\tNA\tNA\t" (q in min_asvs ? min_asvs[q] : "NA");
                                }
                        S   }
                        }
                    ' "${intermediate_dir}/expected_ids.txt" "${intermediate_dir}/recall_blast_all${suffix}.txt" > "${recall_file}"
                else 
                    # NGTax - combine both blast results (forward and reverse) 
                    echo "Combining forward and reverse blast results for NGTax2..."
                    python /projects/0/gs/danique/pipelines/scripts/benchmark/recall_NGTax_asv.py "${intermediate_dir}/recall_blast_all${suffix}.txt" "${intermediate_dir}/recall_blast_all_r${suffix}.txt" > "${recall_file}"

                fi 

                # Get IDs from recall file (excluding header)
                awk 'NR>1{print $1}' "$recall_file" > ${output_dir}/recall_ids${suffix}.txt

                # Add all IDs that are in the expected IDs but not yet in the recall file (i.e. no hits at all)
                awk '
                NR==FNR { recall[$1]=1; next }
                { if (!($1 in recall)) print $1 "\tNA\tNA\tNA" }
                ' ${output_dir}/recall_ids${suffix}.txt "$expected_ids" >> "$recall_file"

                rm -r $DB_DIR

            done
        done
    done
fi
