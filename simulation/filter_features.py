#!/usr/bin/env python3

# Subsets the feature table to the samples kept in the metadata.

import argparse
import pandas as pd

def filter_features(feature_table_file, metadata_file, output_file):
    # Read feature table: rows = features, columns = samples
    features = pd.read_csv(feature_table_file, sep='\t', index_col=0)
    
    # Read metadata: rows = samples
    meta = pd.read_csv(metadata_file, sep='\t', index_col=0, dtype=str)
    samples_to_keep = meta.index.tolist()
    
    # Subset feature table to only keep samples in metadata
    features = features.loc[:, features.columns.intersection(samples_to_keep)]

    # Save filtered features
    features.to_csv(output_file, sep='\t')

def main():
    parser = argparse.ArgumentParser(description="Subset feature table to samples in metadata")
    parser.add_argument('--feature_table', required=True, help='Feature table file (TSV)')
    parser.add_argument('--metadata', required=True, help='Metadata file (TSV)')
    parser.add_argument('--output', required=True, help='Output filtered feature table (TSV)')
    
    args = parser.parse_args()
    
    filter_features(args.feature_table, args.metadata, args.output)

if __name__ == "__main__":
    main()
