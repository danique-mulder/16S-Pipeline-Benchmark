#!/usr/bin/env python3

# Keeps only the samples matching a metadata value (control participants).

import argparse
import pandas as pd

def filter_metadata(metadata_file, output_file, column, value):
    # Read metadata
    meta = pd.read_csv(metadata_file, sep='\t', dtype=str)
    
    if column not in meta.columns:
        raise ValueError(f"Column '{column}' not found in metadata.")
    
    # Filter rows where column == value
    filtered_meta = meta[meta[column] == value]
    
    # Save filtered metadata
    filtered_meta.to_csv(output_file, sep='\t', index=False)

def main():
    parser = argparse.ArgumentParser(description="Filter metadata by column and value.")
    parser.add_argument('--metadata', required=True, help='Input metadata file (TSV)')
    parser.add_argument('--output', required=True, help='Output filtered metadata file (TSV)')
    parser.add_argument('--column', required=True, help='Column to filter on')
    parser.add_argument('--value', required=True, help='Value to filter for in the column')
    
    args = parser.parse_args()
    
    filter_metadata(args.metadata, args.output, args.column, args.value)

if __name__ == "__main__":
    main()
