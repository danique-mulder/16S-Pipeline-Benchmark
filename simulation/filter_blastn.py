#!/usr/bin/env python3
import argparse
import pandas as pd
import matplotlib.pyplot as plt

VAGUE_TERMS = [
    'incertae_sedis', 'incertae sedis', 'unclassified', 'unknown', 'novel',
    'unidentified', 'uncultured', 'metagenome', 'candidate', 'ambiguous',
    'possible', 'probable', 'genomosp', 'phylum', 'class', 'chloroplast', 'mitochondria',
    "ucg-", "gca-", "_group", "family", "vadin", "cag-", "clade"
]

def read_taxonomy(tax_file):
    """Read taxonomy file: seq_id -> taxonomy string"""
    tax = {}
    with open(tax_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                tax[parts[0]] = parts[1]
    return tax

def has_informative_species(tax_str):
    """Return True if taxonomy has a non-vague species label."""
    if not tax_str:
        return False
    for level in tax_str.split(';'):
        if level.strip().startswith("s__"):
            species = level.replace("s__", "").strip().lower()
            if not species:
                return False
            if any(vague in species for vague in VAGUE_TERMS):
                return False
            return True
    return False

def filter_blast(blast_file, taxonomy_file, output_filtered, output_count):
    # Load BLAST results (outfmt 6: qseqid sseqid pident length evalue bitscore)
    cols = ['qseqid','sseqid','pident','length','evalue','bitscore']
    df = pd.read_csv(blast_file, sep="\t", names=cols, header=None)

    # Plot histogram of %identity
    plt.hist(df['pident'], bins=50, color='steelblue', edgecolor='black')
    plt.xlabel("% Identity")
    plt.ylabel("Count")
    plt.title("Histogram of BLAST %identity")
    plt.savefig(output_filtered.replace(".tsv", "_pident_hist.png"))
    plt.close()

    # Filter by identity threshold
    before = len(df)
    df = df[df['pident'] >= 90]
    after = len(df)
    print(f"Filtered by %identity >= 90: kept {after}/{before} hits")

    # Load taxonomy and add as column
    taxonomy = read_taxonomy(taxonomy_file)
    df['taxonomy'] = df['sseqid'].map(taxonomy).fillna("Unknown")

    # Save filtered BLAST results with taxonomy
    blast_file_tax = blast_file.replace(".tsv", "_wtax.tsv")
    df.to_csv(blast_file_tax, sep="\t", index=False)

    # Group by query, select best hit
    best_hits = []
    for qseqid, group in df.groupby('qseqid'):
        # First try species-level informative hits
        species_hits = group[group['taxonomy'].apply(has_informative_species)]
        if not species_hits.empty:
            chosen = species_hits.sort_values('bitscore', ascending=False).iloc[0]
        else:
            chosen = group.sort_values('bitscore', ascending=False).iloc[0]
        best_hits.append(chosen)

    best_df = pd.DataFrame(best_hits)

    # Save filtered BLAST results with taxonomy
    best_df.to_csv(output_filtered, sep="\t", index=False)

    # Count table: ASV, SILVA ID, taxonomy
    count_table = best_df[['qseqid','sseqid','taxonomy']]
    count_table.to_csv(output_count, sep="\t", index=False)

    print(f"Final unique ASVs kept: {count_table['qseqid'].nunique()}")

def main():
    parser = argparse.ArgumentParser(description="Filter BLASTn output to best hit per query")
    parser.add_argument('--blast', required=True, help='BLAST results file (outfmt 6)')
    parser.add_argument('--taxonomy', required=True, help='Taxonomy file (TSV: seq_id \t taxonomy)')
    parser.add_argument('--output', required=True, help='Filtered BLAST results file (TSV, with taxonomy)')
    parser.add_argument('--count', required=True, help='Output count table (ASV, SILVA ID, taxonomy)')
    args = parser.parse_args()
    filter_blast(args.blast, args.taxonomy, args.output, args.count)

if __name__ == "__main__":
    main()
