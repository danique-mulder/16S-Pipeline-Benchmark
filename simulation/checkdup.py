#!/usr/bin/env python3

import argparse
from collections import defaultdict
from Bio import SeqIO
import re

VAGUE_TERMS = [
    'incertae_sedis', 'incertae sedis', 'unclassified', 'unknown', 'novel',
    'unidentified', 'uncultured', 'metagenome', 'candidate', 'ambiguous',
    'possible', 'probable', 'genomosp', 'phylum', 'class', 'chloroplast', 'mitochondria',
    "ucg-", "gca-", "_group", "family", "vadin", "cag-", "clade"
]

def parse_uc(uc_file):
    clusters = defaultdict(list)
    with open(uc_file) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if fields[0] == 'S':
                rep = fields[8]
                clusters[rep].append(rep)
            elif fields[0] == 'H':
                hit = fields[8]
                rep = fields[9]
                clusters[rep].append(hit)
    return clusters

def read_taxonomy(tax_file):
    tax = {}
    with open(tax_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                tax[parts[0]] = parts[1]
    return tax

def extract_genus(tax_str):
    match = re.search(r'g__([^;]+)', tax_str)
    if match:
        return match.group(1).strip()
    return None

def is_vague_taxonomy(taxonomy):
    fields = taxonomy.lower().split(';')[:6]
    return any(term in field for field in fields for term in VAGUE_TERMS)

def has_redundant_labels(taxonomy):
    fields = taxonomy.split(';')
    seen = set()
    for field in fields:
        label = re.sub(r'^[a-z]__', '', field.strip().lower())
        if label in seen:
            return True
        if label:
            seen.add(label)
    return False

def main(args):
    clusters = parse_uc(args.uc)
    taxonomy = read_taxonomy(args.tax)
    sequences = SeqIO.to_dict(SeqIO.parse(args.fasta, "fasta"))

    passed_ids = set()
    filtered_taxonomy = {}
    excluded = []

    total_clusters = len(clusters)
    failed_genus_filter = 0
    failed_vague_taxonomy = 0
    failed_redundant_taxonomy = 0
    failed_seq_ambiguous = 0
    failed_seq_length = 0
    final_passed = 0

    for rep, members in clusters.items():
        genera = set()
        cluster_vague = False
        cluster_redundant = False

        for m in members:
            tax_str = taxonomy.get(m)
            if not tax_str:
                genera.add(None)
                continue

            genus = extract_genus(tax_str)
            genera.add(genus)

            if is_vague_taxonomy(tax_str):
                cluster_vague = True
            if has_redundant_labels(tax_str):
                cluster_redundant = True

        # Exclude cluster if genus inconsistent
        if None in genera or len(genera) != 1:
            failed_genus_filter += 1
            excluded.append((rep, "inconsistent_genus_in_cluster", taxonomy.get(rep, "")))
            continue

        # Exclude cluster if any member vague or redundant
        if cluster_vague:
            failed_vague_taxonomy += 1
            excluded.append((rep, "vague_taxonomy_in_cluster", taxonomy.get(rep, "")))
            continue

        if cluster_redundant:
            failed_redundant_taxonomy += 1
            excluded.append((rep, "redundant_taxonomy_in_cluster", taxonomy.get(rep, "")))
            continue

        # Passed all checks → keep rep
        passed_ids.add(rep)
        filtered_taxonomy[rep] = taxonomy.get(rep)

    # Write outputs
    with open(args.output, 'w') as fasta_out, open(args.tax_out, 'w') as tax_out:
        for seq_id in passed_ids:
            record = sequences.get(seq_id)
            if not record:
                continue
            seq = str(record.seq).upper()

            if any(base not in 'ACGT' for base in seq):
                failed_seq_ambiguous += 1
                excluded.append((seq_id, "ambiguous_bases", filtered_taxonomy[seq_id]))
                continue

            if not (252 <= len(seq) <= 260):
                failed_seq_length += 1
                excluded.append((seq_id, "invalid_length", filtered_taxonomy[seq_id]))
                continue

            SeqIO.write(record, fasta_out, "fasta")
            tax_out.write(f"{seq_id}\t{filtered_taxonomy[seq_id]}\n")
            final_passed += 1

    # Write excluded table
    with open(args.excluded, 'w') as excl_out:
        excl_out.write("sequence_id\treason\ttaxonomy\n")
        for sid, reason, tax in excluded:
            excl_out.write(f"{sid}\t{reason}\t{tax}\n")

    # Print summary
    print("\n=== Filtering Summary ===")
    print(f"Total clusters:                      {total_clusters}")
    print(f"Failed genus consistency:            {failed_genus_filter}")
    print(f"Failed vague taxonomy (any member):  {failed_vague_taxonomy}")
    print(f"Failed redundant taxonomy (any):     {failed_redundant_taxonomy}")
    print(f"Failed due to ambiguous bases:       {failed_seq_ambiguous}")
    print(f"Failed due to length outside 252-260:{failed_seq_length}")
    print(f"Final passed sequences:              {final_passed}")
    print(f"Excluded taxonomy written to:        {args.excluded}")
    print("==============================\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter dereplicated SILVA sequences with unambiguous genus-level labels and remove ambiguous sequences.")
    parser.add_argument("--uc", required=True, help="VSEARCH .uc file from dereplication")
    parser.add_argument("--fasta", required=True, help="Original SILVA FASTA file")
    parser.add_argument("--tax", required=True, help="Taxonomy file (TSV)")
    parser.add_argument("--output", required=True, help="Filtered FASTA output path")
    parser.add_argument("--tax_out", required=True, help="Filtered taxonomy output path")
    parser.add_argument("--excluded", required=True, help="Path to write excluded taxa with reason (TSV format)")
    args = parser.parse_args()
    main(args)
