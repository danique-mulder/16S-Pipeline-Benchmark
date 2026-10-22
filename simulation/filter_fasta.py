#!/usr/bin/env python3
import argparse
from Bio import SeqIO

def filter_fasta(fasta_file, feature_table_file, output_file):
    # Read feature IDs from feature table
    with open(feature_table_file) as f:
        header = f.readline()
        feature_ids = set()
        for line in f:
            feature_ids.add(line.split('\t')[0])
    
    # Filter fasta sequences
    records = (rec for rec in SeqIO.parse(fasta_file, "fasta") if rec.id in feature_ids)
    
    # Write filtered fasta
    count = SeqIO.write(records, output_file, "fasta")
    print(f"Filtered {count} sequences.")

def main():
    parser = argparse.ArgumentParser(description="Filter fasta sequences by feature table IDs")
    parser.add_argument('--fasta', required=True, help='Input fasta file')
    parser.add_argument('--feature_table', required=True, help='Filtered feature table (TSV)')
    parser.add_argument('--output', required=True, help='Output filtered fasta file')
    
    args = parser.parse_args()
    filter_fasta(args.fasta, args.feature_table, args.output)

if __name__ == "__main__":
    main()
