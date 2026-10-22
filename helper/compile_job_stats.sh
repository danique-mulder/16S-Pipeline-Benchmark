#!/bin/bash

# Check if the user provided a log file
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <log_file>"
    exit 1
fi

LOG_FILE="$1"
LOG_DIR="/gpfs/work2/0/gs/danique/pipelines/log"
OUTPUT_FILE="${LOG_DIR}/job_summary_$(basename ${LOG_FILE%.tsv}).tsv"  # Output file name based on input log file

# Check if the log file exists
if [ ! -f "$LOG_FILE" ]; then
    echo "Error: Log file '$LOG_FILE' not found!"
    exit 1
fi

# Header for the output file (excluding date-related fields)
echo -e "Script\tJobID\tWalltime\tElapsed\tCPU_Time\tAve_CPU%\tMax_Mem\tStatus" > "$OUTPUT_FILE"

# Extract job IDs from log file and retrieve stats
awk 'NR>1 {print $2}' "$LOG_FILE" | while read -r jobid; do
    # Get SLURM job stats using sacct
    job_stats=$(sacct -j "$jobid" --format=JobID,Timelimit,Elapsed,TotalCPU,AveCPU,MaxRSS,State%20 --parsable2 --noheader)

    # Clean the job stats by removing any extra newlines and spaces
    job_stats_clean=$(echo "$job_stats" | tr -s '\n' ' ')

    # Extract values (removing any leading/trailing spaces)
    walltime=$(echo "$job_stats_clean" | awk -F '|' '{print $2}' | xargs)
    elapsed=$(echo "$job_stats_clean" | awk -F '|' '{print $3}' | xargs)
    cpu_time=$(echo "$job_stats_clean" | awk -F '|' '{print $4}' | xargs)
    ave_cpu=$(echo "$job_stats_clean" | awk -F '|' '{print $5}' | xargs)
    max_mem=$(echo "$job_stats_clean" | awk -F '|' '{print $6}' | xargs)
    status=$(echo "$job_stats_clean" | awk -F '|' '{print $7}' | xargs)

    # Get script name from log file
    script_name=$(awk -v jobid="$jobid" '$2 == jobid {print $1}' "$LOG_FILE")

    # Append to output file, ensuring proper formatting
    echo -e "${script_name}\t${jobid}\t${walltime}\t${elapsed}\t${cpu_time}\t${ave_cpu}\t${max_mem}\t${status}" >> "$OUTPUT_FILE"
    
done

cat "$OUTPUT_FILE"
echo "Job stats compiled into: $OUTPUT_FILE"
