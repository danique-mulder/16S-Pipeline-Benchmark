#!/usr/bin/env python3

#  Sums the counts of ASVs that map to the same representative sequence to build the simulated count table.

import argparse
import pandas as pd
from collections import defaultdict

def main(args):
    # Read count table (ASV to SILVA mapping + taxonomy)
    count_df = pd.read_csv(args.count, sep='\t')
    # Expected columns: sseqid, qseqid, taxonomy

    # Build mappings
    asv_to_silva = dict(zip(count_df['qseqid'], count_df['sseqid']))
    silva_to_tax = dict(zip(count_df['sseqid'], count_df['taxonomy']))

    # Group ASVs per SILVA ID for mapping output
    silva_to_asvs = defaultdict(list)
    for asv, silva_id in asv_to_silva.items():
        silva_to_asvs[silva_id].append(asv)

    # Read feature table (ASVs as rows, samples as columns)
    feature_df = pd.read_csv(args.feature_table, sep='\t', index_col=0)

    # Map feature table rows (ASVs) to SILVA IDs
    # Filter to ASVs present in the count table mapping
    mapped_asvs = [asv for asv in feature_df.index if asv in asv_to_silva]
    feature_df_filtered = feature_df.loc[mapped_asvs]

    # Create a new dataframe for merged counts indexed by SILVA IDs
    merged_df = pd.DataFrame(0, index=pd.Index(list(silva_to_asvs.keys()), name='sseqid'), columns=feature_df.columns)

    # Aggregate counts: sum counts for all ASVs that map to the same SILVA ID
    for silva_id, asvs in silva_to_asvs.items():
        # sum counts across ASVs belonging to this silva_id
        subset = feature_df_filtered.loc[asvs]
        # If only one ASV, subset will be Series, convert to DataFrame
        if isinstance(subset, pd.Series):
            subset = subset.to_frame().T
        merged_df.loc[silva_id] = subset.sum(axis=0)

    # Rename the index column to SILVA_ID
    # merged_df.index.name = 'SILVA_ID'

    # Write merged counts table
    merged_df.to_csv(args.output, sep='\t')

    # Write mapping file: SILVA ID -> semicolon-separated ASVs
    with open(args.mapping, 'w') as f:
        f.write("SILVA_ID\tASV_ID\n")
        for silva_id, asvs in silva_to_asvs.items():
            f.write(f"{silva_id}\t{';'.join(asvs)}\n")

    # Write taxonomy file: SILVA ID -> taxonomy
    with open(args.taxonomy, 'w') as f:
        f.write("SILVA_ID\tTaxon\n")
        for silva_id, tax in silva_to_tax.items():
            f.write(f"{silva_id}\t{tax}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Merge count table with feature table mapping ASVs to SILVA IDs.")
    parser.add_argument("--count", required=True, help="Count table with columns: sseqid, qseqid, taxonomy")
    parser.add_argument("--feature_table", required=True, help="Feature table with ASVs as rows and sample counts as columns")
    parser.add_argument("--mapping", required=True, help="Output path for mapping file (SILVA ID to ASVs)")
    parser.add_argument("--taxonomy", required=True, help="Output path for taxonomy file (SILVA ID to taxonomy)")
    parser.add_argument("--output", required=True, help="Output path for merged counts table")
    args = parser.parse_args()

    main(args)
