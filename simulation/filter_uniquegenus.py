#!/usr/bin/env python3

#  Keeps only representatives with a single, unambiguous genus.

import argparse
import pandas as pd
from Bio import SeqIO

def main(args):
    # Step 1: Read alt_classifications and filter for unique genus
    alt_class = pd.read_csv(args.alt_class, sep="\t")
    resolvable_ids = set(alt_class.loc[alt_class["num_genera"] == 1, "seq_id"])

    print(f"Total ASVs in alt_class: {len(alt_class)}")
    print(f"Resolvable ASVs (unique genus): {len(resolvable_ids)}")

    # Step 2: Subset count_table
    count_table = pd.read_csv(args.count_table, sep="\t", index_col=1)
    filtered_count = count_table.loc[count_table.index.intersection(resolvable_ids)]
    filtered_count.to_csv(args.output_count, sep="\t")

    print(f"Count table filtered: {filtered_count.shape[0]} ASVs retained")

    # Step 3: Subset fasta
    with open(args.output_fasta, "w") as out_fasta:
        retained = 0
        for record in SeqIO.parse(args.tosimulate, "fasta"):
            if record.id in resolvable_ids:
                SeqIO.write(record, out_fasta, "fasta")
                retained += 1

    print(f"FASTA filtered: {retained} sequences written")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Keep only ASVs with unique genus classification"
    )
    parser.add_argument("--count_table", required=True, help="Count table TSV")
    parser.add_argument("--alt_class", required=True, help="Alternative classifications TSV")
    parser.add_argument("--tosimulate", required=True, help="Input FASTA with sequences")
    parser.add_argument("--output_count", required=True, help="Output filtered count table")
    parser.add_argument("--output_fasta", required=True, help="Output filtered FASTA")
    args = parser.parse_args()
    main(args)
