#!/bin/bash
#SBATCH -p rome
#SBATCH -t 5:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output=/gpfs/work2/0/gs/danique/pipelines/data/remux_output_%j.out

# Re-multiplexes demultiplexed reads into one R1/R2 file using artificial 8-nt barcodes, and writes the barcode file. 
# Used for eDNAFlow.

while getopts ":i:p:o:f:r:h:" opt; do
  case $opt in
    i) BASE_DIR="$OPTARG";;
    p) PRIMER_DIR="$OPTARG";;
    o) BARCODE_DIR="$OPTARG";;
    f) FORWARD_PRIMER="$OPTARG";;
    r) REVERSE_PRIMER="$OPTARG";;
    h)
      echo "Usage: \$0 [options]"
      echo "Options:"
      echo "  -i <base_dir>           Specify base directory."
      echo "  -p <primer_dir>         Specify primer directory. This directory should contain all fastq.gz with primers attached"
      echo "  -f <forward_primer>     Specify forward primer sequence"
      echo "  -r <reverse_primer>     Specify reverse primer sequence"
      echo "  -h                      Display this help message"
      exit 0
      ;;
    \?) echo "Invalid option: -$OPTARG"; exit 1;;
  esac
done

cd $BASE_DIR

# New directory
if ! [ -d "${BASE_DIR}/remux" ]; then
  mkdir "${BASE_DIR}/remux"
fi

echo "Creating barcodes..."

# Create the header for the barcodes.tsv file
echo -e "#ASSAY\tSAMPLE_ID\tF_TAG_SEQUENCE:R_TAG_SEQUENCE\tFORWARD_PRIMER\tREVERSE_PRIMER" > barcodes.txt

# Declare an associative array to store unique barcodes
declare -A barcode_set

# Loop through the files and generate unique barcodes
for pathname in $PRIMER_DIR/*R1*.gz; do
    filename=$(basename "$pathname") # extract filename
    sample=$(echo "$filename" | awk -F'_' '{print $1}') # extract sample name

    echo "Processing ${sample}..."

    # Generate a unique 8-nt barcode
    while true; do
        barcode=$(printf "%08d" $((RANDOM % 100000000)) | tr '0123456789' 'ACGT')
        if [[ -z "${barcode_set[$barcode]}" ]]; then
            barcode_set["$barcode"]=1  # Mark this barcode as used
            break
        fi
    done

    # Append the assay, sample name, and barcode to the file
    echo -e "16S\t${sample}\t${barcode}:${barcode}\t${FORWARD_PRIMER}\t${REVERSE_PRIMER}" >> barcodes.txt
done

echo
head barcodes.txt
echo

# Initialize the output files (R1 and R2)
f_file="${BASE_DIR}/seq/remux/remux_R1.fastq"
> "$f_file"  # Clear the file if it already exists

r_file="${BASE_DIR}/seq/remux/remux_R2.fastq"
> "$r_file"  # Clear the file if it already exists

# Loop through the forward reads and append the barcodes
for pathname_f in $PRIMER_DIR/*R1*.gz; do
    filename=$(basename "$pathname_f") # extract filename
    samplename=$(echo "$filename" | awk -F'_' '{print $1}') # extract sample name
    
    echo "processing ${filename}..."
    
    # retrieve barcode
    barcode=$(awk -v sample="$samplename" '$2 == sample {print $3}' barcodes.txt)
    barcode=$(echo $barcode | awk -F':' '{print $1}')
    
    barcode_quality=$(printf 'I%.0s' $(seq 1 ${#barcode})) # generate (high) quality scores
    
    zcat "$pathname_f" | awk -v barcode="$barcode" -v barcode_quality="$barcode_quality" '
     {
       if (NR % 4 == 1) {
           print $0 # header
       } else if (NR % 4 == 2) {
           print barcode $0
       } else if (NR % 4 == 3) {
           print $0
       } else if (NR % 4 == 0) {
           print barcode_quality $0
       }
     }' >> "$f_file"
    
  pathname_r="${pathname_f/R1/R2}"
  if [ ! -f "$pathname_r" ]; then
    echo "Reverse read file not found for $pathname_f" >&2
  else
    zcat "$pathname_r" | awk -v barcode="$barcode" -v barcode_quality="$barcode_quality" '
     {
          if (NR % 4 == 1) {
              print $0 # header
          } else if (NR % 4 == 2) {
              print barcode $0
          } else if (NR % 4 == 3) {
              print $0
          } else if (NR % 4 == 0) {
              print barcode_quality $0
          }
    }' >> "$r_file"
  fi

done

# Initialize a variable to store the total line count
total_lines=0

# Loop through all R1 files and count lines
for pathname in ${PRIMER_DIR}/*R1*.gz; do
  # Count the number of lines in the current file
  file_lines=$(zcat "$pathname" | wc -l)
  
  # Add the count to the total
  total_lines=$((total_lines + file_lines))
  
  # Print the line count for the current file (optional)
  echo "File: $pathname          Lines: $file_lines"
done

# Print the total number of lines across all files
echo "Total lines across all files (forward):  $total_lines"

remux_lines=$(cat "$f_file" | wc -l)
echo "Total lines in remux file (forward):     $remux_lines"


# Initialize a variable to store the total line count
total_lines=0

# Loop through all R1 files and count lines
for pathname in ${PRIMER_DIR}/*R2*.gz; do
  # Count the number of lines in the current file
  file_lines=$(zcat "$pathname" | wc -l)
  
  # Add the count to the total
  total_lines=$((total_lines + file_lines))
  
  # Print the line count for the current file (optional)
  echo "File: $pathname,          Lines: $file_lines"
done

# Print the total number of lines across all files
echo "Total lines across all files (forward):  $total_lines"

remux_lines=$(cat "$r_file" | wc -l)
echo "Total lines in remux file (reverse):     $remux_lines"
