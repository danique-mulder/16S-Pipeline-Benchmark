#!/bin/bash

SCRIPT_DIR="/gpfs/work2/0/gs/danique/pipelines/scripts/pipelines"
LOG_DIR="/gpfs/work2/0/gs/danique/pipelines/log"
cd "$SCRIPT_DIR"

# Function to run pipelines
run_pipelines() {

    # Reset OPTIND to ensure getopts works for each function call
    OPTIND=1

    # Parse command-line options
    while getopts "d:i:b:a:f:r:p:m:l:s:c:t:u:e:j:k:z:x:y:n:g:h" opt; do
        case $opt in
            d) DATASET="$OPTARG" ;;
            i) DATA_DIR="$OPTARG" ;;
            b) BARCODEFILE="$OPTARG" ;;
            a) PRIMER_ATTACHED="$OPTARG" ;;
            f) FORWARD_PRIMER="$OPTARG" ;;
            r) REVERSE_PRIMER="$OPTARG" ;;
            p) REVCOMP="$OPTARG" ;;
            m) MACHINE="$OPTARG" ;;
            l) SEQUENCING_RUN="$OPTARG" ;;
            s) SEPARATOR="$OPTARG" ;;
            c) CLUSTER="$OPTARG" ;;
            t) TRUNC_F="$OPTARG" ;;
            u) TRUNC_R="$OPTARG" ;;
            e) PATTERN="$OPTARG" ;;
            j) MINOVLEN="$OPTARG" ;;
            k) MAXDIFF="$OPTARG" ;;
            z) RUN_SHORT="$OPTARG" ;;
            x) RUN_LONG="$OPTARG" ;;
            y) REMOVEOLD="$OPTARG" ;;
            n) TASKS="$OPTARG" ;;
            g) THREADS_PER_TASK="$OPTARG" ;;
            h) echo "Usage: $0 [-d DATASET] [-i DATA_DIR] [-b BARCODEFILE] [-a PRIMER_ATTACHED] [-f FORWARD_PRIMER] [-r REVERSE_PRIMER] [-p REVCOMP] [-m MACHINE] [-l SEQUENCING_RUN] [-s SEPARATOR] [-c CLUSTER] [-t TRUNC_F] [-u TRUNC_R] [-e PATTERN] [-z RUN_SHORT] [-x RUN_LONG] [-y REMOVEOLD] [-n TASKS] [-g THREADS_PER_TASK]" >&2; exit 1 ;;
            *) echo "Invalid option: -$opt" >&2; exit 1 ;;
        esac
    done

    if [ $PRIMER_ATTACHED -eq 0 ]; then
        local DATA_DIR_PRIMER="${DATA_DIR}/primer"
    else
        local DATA_DIR_PRIMER="${DATA_DIR}"
    fi

    local DATA_DIR_REMUX="${DATA_DIR}/remux"
    # Use timestamp to avoid overwriting log files for the same dataset
    local LOG_FILE="${LOG_DIR}/log_${DATASET}_$(date +%Y%m%d_%H%M%S).tsv"
    local THREADS=$((TASKS * THREADS_PER_TASK))

    declare -A pipelines=(
        ["AMPtk.sh"]="-i ${DATA_DIR} -a ${PRIMER_ATTACHED} -y ${REMOVEOLD}"
        ["VSEARCH.sh"]="-i ${DATA_DIR} -a ${PRIMER_ATTACHED} -j ${MINOVLEN} -k ${MAXDIFF} -y ${REMOVEOLD}"
        ["DADA2.sh"]="-i ${DATA_DIR} -a ${PRIMER_ATTACHED} -t ${TRUNC_F} -u ${TRUNC_R}"
        ["Tourmaline.sh"]="-i ${DATA_DIR} -a ${PRIMER_ATTACHED} -t ${TRUNC_F} -u ${TRUNC_R} -e ${PATTERN}"
        ["dadasnake.sh"]="-i ${DATA_DIR} -a ${PRIMER_ATTACHED} -m ${MACHINE} -t ${TRUNC_F} -u ${TRUNC_R} -y ${REMOVEOLD}"
        ["NGTax2.sh"]="-i ${DATA_DIR} -a ${PRIMER_ATTACHED} -s ${SEPARATOR} -y ${REMOVEOLD}"
        ["nfcore_ampliseq.sh"]="-i ${DATA_DIR} -a ${PRIMER_ATTACHED} -s ${SEPARATOR} -y ${REMOVEOLD}"
        ["Natrix.sh"]="-i ${DATA_DIR} -a ${PRIMER_ATTACHED} -c asv -y ${REMOVEOLD} -g ${THREADS}" 
        ["Natrix.sh"]="-i ${DATA_DIR} -a ${PRIMER_ATTACHED} -c otu -y ${REMOVEOLD} -g ${THREADS}"
        ["DADAist2.sh"]="-i ${DATA_DIR} -a ${PRIMER_ATTACHED} -t ${TRUNC_F} -u ${TRUNC_R} -y ${REMOVEOLD}"
        ["QIIME2.sh"]="-i ${DATA_DIR} -a ${PRIMER_ATTACHED} -t ${TRUNC_F} -u ${TRUNC_R} -y ${REMOVEOLD}"
        ["QIIME2_deblur.sh"]="-i ${DATA_DIR} -a ${PRIMER_ATTACHED} -y ${REMOVEOLD}"
        ["MetaWorks.sh"]="-i ${DATA_DIR_PRIMER} -p ${REVCOMP} -a 1 -t ${TRUNC_F} -u ${TRUNC_R} -y ${REMOVEOLD}"
        ["APSCALE.sh"]="-i ${DATA_DIR_PRIMER} -a 1 -y ${REMOVEOLD}"
        ["LotuS2.sh"]="-i ${DATA_DIR_PRIMER} -a 1 -m ${MACHINE} -l ${SEQUENCING_RUN} -c asv" 
        ["LotuS2.sh"]="-i ${DATA_DIR_PRIMER} -a 1 -m ${MACHINE} -l ${SEQUENCING_RUN} -c otu"
        ["eDNAFlow.sh"]="-i ${DATA_DIR_REMUX} -b ${BARCODEFILE} -y 0 -g ${THREADS}"
    )

    # Create a new log file
    echo -e "#script\tjobid\tcommand" > "${LOG_FILE}"

    # Loop over the pipelines
    for script in "${!pipelines[@]}"; do
        # Assign runtime based on the script
        if [[ "$script" == "OBITOOLS.sh" || "$script" == "eDNAFlow.sh" || "$script" == "AMPtk.sh" ]]; then
            RUNTIME=$RUN_LONG
        else
            RUNTIME=$RUN_SHORT
        fi

        # If multiple configurations exist, split them by the pipe character
        IFS='|' read -r -a options <<< "${pipelines[$script]}"

        # Loop over the configurations for the script
        for option in "${options[@]}"; do
            command="sbatch -p genoa -t ${RUNTIME} -n ${TASKS} -c ${THREADS_PER_TASK} $script -d ${DATASET} -f ${FORWARD_PRIMER} -r ${REVERSE_PRIMER} $option"
    
            # Execute sbatch and capture job ID
            jobid_output=$($command 2>&1)
            if [[ $jobid_output =~ Submitted\ batch\ job\ ([0-9]+) ]]; then
                jobid=${BASH_REMATCH[1]}
            else
                jobid="FAILED"
            fi
            # Log the command and jobid
            echo -e "${script}\t${jobid}\t${command}" >> "${LOG_FILE}"
        done
    done

    cat "${LOG_FILE}"
}

###################################################################################################################################################
#                                                                 Running pipelines                                                               #
###################################################################################################################################################

dats="parkinson" # mock5 mock12 gemsim mindset parkinson miseq_v4_sim novaseq_v4_sim mock5_o mock12_o

# Run pipelines for mock5
if [[ $dats == *mock5* ]]; then
    run_pipelines \
        -d "mock5" \
        -i "/gpfs/work2/0/gs/danique/mockrobiota/mock5/fastq_demux" \
        -b "/gpfs/work2/0/gs/danique/mockrobiota/mock5/fastq_demux/barcodes.txt" \
        -a 0 \
        -f "GTGCCAGCMGCCGCGGTAA" \
        -r "GGACTACHVGGGTWTCTAAT" \
        -p "ATTAGAWACCCBDGTAGTCC" \
        -m "miSeq" \
        -l 0 \
        -s "_" \
        -c "asv" \
        -t 190 \
        -u 135 \
        -e "_L001_R[12]_001.fastq.gz" \
        -j 200 \
        -k 15 \
        -z "2:00:00" \
        -x "12:00:00" \
        -y 1 \
        -n 1 \
        -g 24
fi

# Run pipelines for mock12
if [[ $dats == *mock12* ]]; then
    run_pipelines \
        -d "mock12" \
        -i "/gpfs/work2/0/gs/danique/mockrobiota/mock12/fastq_demux" \
        -b "/gpfs/work2/0/gs/danique/mockrobiota/mock12/fastq_demux/barcodes.txt" \
        -a 0 \
        -f "GTGCCAGCMGCCGCGGTAA" \
        -r "GGACTACHVGGGTWTCTAAT" \
        -p "ATTAGAWACCCBDGTAGTCC" \
        -m "miSeq" \
        -l 0 \
        -s "_" \
        -c "asv" \
        -t 190 \
        -u 135 \
        -e "_L001_R[12]_001.fastq.gz" \
        -j 200 \
        -k 15 \
        -z "2:00:00" \
        -x "12:00:00" \
        -y 1 \
        -n 1 \
        -g 24
fi

# Run pipelines for gemsim
if [[ $dats == *gemsim* ]]; then
    run_pipelines \
        -d "gemsim" \
        -i "/gpfs/work2/0/gs/danique/pipelines/data/gemsim/seq" \
        -a 1 \
        -f "GTGCCAGCMGCCGCGGTAA" \
        -r "GGACTACHVGGGTWTCTAAT" \
        -p "ATTAGAWACCCBDGTAGTCC" \
        -m "novaSeq" \
        -l 0 \
        -s "_" \
        -c "asv" \
        -t 230 \
        -u 230 \
        -e "_L001_R[12]_001.fastq.gz" \
        -j 200 \
        -k 15 \
        -z "5-00:00:00" \
        -x "5-00:00:00" \
        -y 1 \
        -n 1 \
        -g 48
fi

# Run pipelines for mindset
if [[ $dats == *mindset* ]]; then
    run_pipelines \
        -d "mindset" \
        -i "/gpfs/work2/0/gs/danique/pipelines/data/mindset/seq" \
        -a 1 \
        -f "GTGCCAGCMGCCGCGGTAA" \
        -r "GGACTACHVGGGTWTCTAAT" \
        -p "ATTAGAWACCCBDGTAGTCC" \
        -m "novaSeq" \
        -l 0 \
        -s "_" \
        -c "asv" \
        -t 230 \
        -u 230 \
        -e "_L001_R[12]_001.fastq.gz" \
        -j 200 \
        -k 15 \
        -z "2-00:00:00" \
        -x "5-00:00:00" \
        -y 1 \
        -n 1 \
        -g 48
fi

# Run pipelines for parkinson
if [[ $dats == *parkinson* ]]; then  
    run_pipelines \
        -d "parkinson" \
        -i "/gpfs/work2/0/gs/danique/pipelines/data/parkinson/seq" \
        -a 1 \
        -f "GTGCCAGCMGCCGCGGTAA" \
        -r "GGACTACHVGGGTWTCTAAT" \
        -p "ATTAGAWACCCBDGTAGTCC" \
        -m "miSeq" \
        -l 0 \
        -s "_" \
        -c "asv" \
        -t 230 \
        -u 215 \
        -e "_L001_R[12]_001.fastq.gz" \
        -j 200 \
        -k 15 \
        -z "2-00:00:00" \
        -x "5-00:00:00" \
        -y 1 \
        -n 1 \
        -g 48
fi

# Run pipelines for miseq_v4_sim
if [[ $dats == *miseq_v4_sim* ]]; then    
    run_pipelines \
        -d "miseq_v4_sim" \
        -i "/gpfs/work2/0/gs/danique/pipelines/simulations/input/fastq/miseq" \
        -a 0 \
        -f "GTGCCAGCMGCCGCGGTAA" \
        -r "GGACTACHVGGGTWTCTAAT" \
        -p "ATTAGAWACCCBDGTAGTCC" \
        -m "miSeq" \
        -l 0 \
        -s "_" \
        -c "asv" \
        -t 195 \
        -u 150 \
        -e "_L001_R[12]_001.fastq.gz" \
        -j 200 \
        -k 15 \
        -z "5-00:00:00" \
        -x "5-00:00:00" \
        -y 1 \
        -n 1 \
        -g 48
fi


# Run pipelines for novaseq_v4_sim
if [[ $dats == *novaseq_v4_sim* ]]; then   
    run_pipelines \
        -d "novaseq_v4_sim" \
        -i "/gpfs/work2/0/gs/danique/pipelines/simulations/input/fastq/novaseq" \
        -a 0 \
        -f "GTGCCAGCMGCCGCGGTAA" \
        -r "GGACTACHVGGGTWTCTAAT" \
        -p "ATTAGAWACCCBDGTAGTCC" \
        -m "novaSeq" \
        -l 0 \
        -s "_" \
        -c "asv" \
        -t 150 \
        -u 150 \
        -e "_L001_R[12]_001.fastq.gz" \
        -j 40 \
        -k 3 \
        -z "5-00:00:00" \
        -x "5-00:00:00" \
        -y 1 \
        -n 1 \
        -g 48
fi

### Only for running Natrix with original setting 
# # Run pipelines for mock5
# if [[ $dats == *mock5_o* ]]; then  
#     run_pipelines \
#         -d "mock5_o" \
#         -i "/gpfs/work2/0/gs/danique/mockrobiota/mock5/fastq_demux" \
#         -b "/gpfs/work2/0/gs/danique/mockrobiota/mock5/fastq_demux/barcodes.txt" \
#         -a 0 \
#         -f "GTGCCAGCMGCCGCGGTAA" \
#         -r "GGACTACHVGGGTWTCTAAT" \
#         -p "ATTAGAWACCCBDGTAGTCC" \
#         -m "miSeq" \
#         -l 0 \
#         -s "_" \
#         -c "asv" \
#         -t 190 \
#         -u 135 \
#         -e "_L001_R[12]_001.fastq.gz" \
#         -j 200 \
#         -k 15 \
#         -z "2:00:00" \
#         -x "12:00:00" \
#         -y 1 \
#         -n 1 \
#         -g 24
# fi

# # Run pipelines for mock12
# if [[ $dats == *mock12_o* ]]; then  
#     run_pipelines \
#         -d "mock12_o" \
#         -i "/gpfs/work2/0/gs/danique/mockrobiota/mock12/fastq_demux" \
#         -b "/gpfs/work2/0/gs/danique/mockrobiota/mock12/fastq_demux/barcodes.txt" \
#         -a 0 \
#         -f "GTGCCAGCMGCCGCGGTAA" \
#         -r "GGACTACHVGGGTWTCTAAT" \
#         -p "ATTAGAWACCCBDGTAGTCC" \
#         -m "miSeq" \
#         -l 0 \
#         -s "_" \
#         -c "asv" \
#         -t 190 \
#         -u 135 \
#         -e "_L001_R[12]_001.fastq.gz" \
#         -j 200 \
#         -k 15 \
#         -z "2:00:00" \
#         -x "12:00:00" \
#         -y 1 \
#         -n 1 \
#         -g 24
# fi
