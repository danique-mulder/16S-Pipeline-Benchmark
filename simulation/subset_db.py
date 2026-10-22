#!/usr/bin/env python3

import argparse
from Bio import SeqIO

def main(args):
    # Read list of SILVA IDs to keep from count_table.tsv
    silva_ids = set()
    with open(args.count) as f:
        for line in f:
            if line.startswith("qseqid"):  # skip header if present
                continue
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                silva_id = parts[1]
                silva_ids.add(silva_id)

    # Read SILVA fasta and write only sequences in silva_ids
    with open(args.output, 'w') as out_fasta:
        for record in SeqIO.parse(args.db, 'fasta'):
            if record.id in silva_ids:
                SeqIO.write(record, out_fasta, 'fasta')

    print(f"Subset database written to {args.output} with {len(silva_ids)} sequences.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Subset SILVA fasta database to sequences in filtered blast results")
    parser.add_argument('--db', required=True, help="Filtered SILVA fasta database")
    parser.add_argument('--count', required=True, help="Count table with SILVA IDs (TSV from filter_blastn.py)")
    parser.add_argument('--output', required=True, help="Output fasta with subset sequences")
    args = parser.parse_args()

    main(args)
