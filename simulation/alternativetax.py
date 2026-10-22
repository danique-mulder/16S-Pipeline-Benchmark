#!/usr/bin/env python3
import argparse
import re
from collections import defaultdict

def read_taxonomy(tax_file):
    """Read taxonomy file into dict: seq_id -> taxonomy string."""
    tax = {}
    with open(tax_file) as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) == 2:
                tax[parts[0]] = parts[1]
    return tax

def extract_genus(tax_str):
    """Extract genus label from taxonomy string """
    match = re.search(r"g__([^;]+)", tax_str)
    if match:
        return match.group(1).strip()
    return None

def main(args):
    taxonomy = read_taxonomy(args.taxonomy)
    results = defaultdict(set)

    # Parse BLAST results
    with open(args.blast) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 2:
                continue

            query_id, subject_id = parts[0], parts[1]
            tax_str = taxonomy.get(subject_id, "")
            genus = extract_genus(tax_str)
            if genus:
                results[query_id].add(genus)

    # Write output
    with open(args.output, "w") as out:
        out.write("seq_id\tnum_genera\tgenera\n")
        for query, genera in results.items():
            genera_list = sorted(genera)
            out.write(f"{query}\t{len(genera_list)}\t{';'.join(genera_list)}\n")

    print(f"Processed {len(results)} queries")
    print(f"Output written to {args.output}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract all possible genus-level assignments per query from BLASTn results"
    )
    parser.add_argument("--blast", required=True, help="BLASTn output table (outfmt 6, TSV)")
    parser.add_argument("--taxonomy", required=True, help="Taxonomy file (TSV: ID \\t taxonomy)")
    parser.add_argument("--output", required=True, help="Output TSV with genus-level assignments")
    args = parser.parse_args()
    main(args)
