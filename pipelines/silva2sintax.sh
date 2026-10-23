#!/bin/bash

# Converts a SILVA database to SINTAX format for use with VSEARCH.

# Set variables
# Set default values

while getopts ":d:h:" opt; do
  case $opt in
    d) DB="$OPTARG";;
    h)
      echo "Usage: \$0 [options]"
      echo "Options:"
      echo "  -d <database>    Specify path to database"
      echo "  -o <outputpath>  Specify data directory. This directory should contain all fastq.gz files to be processed"
      echo "  -h               Display this help message"
      exit 0
      ;;
    \?) echo "Invalid option: -$OPTARG"; exit 1;;
  esac
done

OUTDIR=$(dirname "${DB}")
DB_name=$(basename "${DB}")
REFORM="${OUTDIR}/${DB_name}_reformatted.fasta"
SINTAX="${OUTDIR}/${DB_name}_sintax.fasta"

# Reformat so that all odd lines are sequence identifiers (starting with >) and all even lines contain the sequence and so that all occurances of U (uracil) are replaced with T (thyamine)
zcat "${DB}" | \
awk '/^>/{if (NR > 1) print seq; print; seq=""; next} {seq=seq""$0} END{print seq}' | \
awk 'NR % 2 == 0 {gsub(/U/, "T")} {print}' > "${REFORM}"

# Reformat the sequence headers to SINTAX format
cat "${REFORM}" | \
awk '
BEGIN { FS=";"; OFS="," }
# Process header lines
/^>/ {
    # Extract the sequence ID and taxonomy
    split($1, id, " ")
    # Start building the new header
    new_header = id[1] ";tax="
    # Process taxonomy fields
    for (i = 2; i <= NF; i++) {
        # Map taxonomy levels to SINTAX format
        if (i == 2) new_header = new_header "d:" $i
        else if (i == 3 && $i != "") new_header = new_header ",p:" $i
        else if (i == 4 && $i != "") new_header = new_header ",c:" $i
        else if (i == 5 && $i != "") new_header = new_header ",o:" $i
        else if (i == 6 && $i != "") new_header = new_header ",f:" $i
        else if (i == 7 && $i != "") new_header = new_header ",g:" $i
        else if (i == 8 && $i != "") new_header = new_header ",s:" $i
    }
    print new_header ";"
    next
}
# Print sequence lines as they are
{ print }
' > "${SINTAX}"

rm "${REFORM}"
