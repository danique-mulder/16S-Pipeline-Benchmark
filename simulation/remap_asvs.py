
# Map ASVs that were lost during prior steps to the best fitting representative sequence (SILVA ID) based on BLAST results and taxonomy

import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("--blast", required=True, help="BLAST results")
parser.add_argument("--taxonomy", required=True, help="Silva taxonomy file (tsv with ID and taxonomy)")
parser.add_argument("--count", required=True, help="count table (resolvable)")
parser.add_argument("--output", required=True, help="Output remapped table")
args = parser.parse_args()

# Load BLAST results
cols = ["qseqid", "sseqid", "pident", "length", "evalue", "bitscore"]
blast = pd.read_csv(args.blast, sep="\t", names=cols)

# Keep only hits with ≥95% identity
blast = blast[blast["pident"] >= 90]

# Keep only best hit per query (highest bitscore, then highest pident)
blast = blast.sort_values(["qseqid", "bitscore", "pident"], ascending=[True, False, False])
blast_best = blast.drop_duplicates(subset="qseqid", keep="first")

# Load taxonomy
tax = pd.read_csv(args.taxonomy, sep="\t", header=None, names=["sseqid", "taxonomy"])

# Merge blast hits with taxonomy
merged = blast_best.merge(tax, on="sseqid", how="left")

# Load the existing count table
count = pd.read_csv(args.count, sep="\t")

# Concatenate count table + new rows
final = pd.concat([count, merged], ignore_index=True)

# Save output: sseqid, qseqid, taxonomy
final[["sseqid", "qseqid", "taxonomy"]].to_csv(args.output, sep="\t", index=False)