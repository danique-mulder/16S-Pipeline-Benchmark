#!/usr/bin/env python3
"""
precision_NGTax2_asv.py

Combines forward and reverse BLAST results to classify ASVs by total mismatches.

Usage:
    python precision_NGTax2_asv.py forward_blast.tsv reverse_blast.tsv intermediate_dir

Outputs:
    - {intermediate_dir}/smallblast_good.tsv
      ASVs where forward and reverse map to the same reference with total mismatches 0 or 1.

    - {intermediate_dir}/smallblast_bad.tsv
      ASVs where total mismatches ≥ 2.

Assumptions:
- BLAST TSV columns:
  qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
- Uses mismatch from column 5 directly (=mismatch).
- Considers only best (lowest mismatch) hit per ASV per strand.
- Only pairs hits where forward and reverse map to the same sseqid.
"""

import sys
import os
import csv
from collections import defaultdict

def parse_blast_hits(blastfile):
    """
    Parse BLAST TSV, return dict: hits[qid] = dict { sseqid: min_mismatches }
    Uses mismatch from column 5 (0-based index 4).
    """
    hits = defaultdict(lambda: defaultdict(lambda: float('inf')))
    with open(blastfile) as f:
        reader = csv.reader(f, delimiter='\t')
        for row in reader:
            if len(row) < 12:
                continue
            qseqid, sseqid = row[0], row[1]
            mismatches = int(row[4])  # mismatch column

            if mismatches < hits[qseqid][sseqid]:
                hits[qseqid][sseqid] = mismatches
    return hits

def main():
    if len(sys.argv) != 4:
        print("Usage: python precision_NGTax2_asv.py forward.tsv reverse.tsv intermediate_dir")
        sys.exit(1)

    fwd_file = sys.argv[1]
    rev_file = sys.argv[2]
    intermediate_dir = sys.argv[3]

    if not os.path.exists(intermediate_dir):
        os.makedirs(intermediate_dir)

    fwd_hits = parse_blast_hits(fwd_file)
    rev_hits = parse_blast_hits(rev_file)

    good_path = os.path.join(intermediate_dir, "smallblast_good.tsv")
    bad_path = os.path.join(intermediate_dir, "smallblast_bad.tsv")

    with open(good_path, "w", newline='') as good_fh, open(bad_path, "w", newline='') as bad_fh:
        good_writer = csv.writer(good_fh, delimiter='\t')
        bad_writer = csv.writer(bad_fh, delimiter='\t')

        good_writer.writerow(["ID", "tax", "mismatch", "reference"])
        # bad_writer.writerow(["ID", "tax", "mismatch", "reference"])

        all_asvs = set(fwd_hits.keys()) | set(rev_hits.keys())

        for asv in sorted(all_asvs):
            fwd_asv_hits = fwd_hits.get(asv, {})
            rev_asv_hits = rev_hits.get(asv, {})

            # Common sseqids between forward and reverse hits
            common_sseqids = set(fwd_asv_hits.keys()) & set(rev_asv_hits.keys())
            if not common_sseqids:
                # No shared hit in both strands; skip or handle separately if needed
                continue

            combined = []
            for sseqid in common_sseqids:
                total_mismatch = fwd_asv_hits[sseqid] + rev_asv_hits[sseqid]
                combined.append((sseqid, total_mismatch))

            # Select best (lowest mismatch) combined hit for ASV
            best_sseqid, best_mismatch = min(combined, key=lambda x: x[1])

            if best_mismatch <= 1:
                good_writer.writerow([asv, best_sseqid, best_mismatch, "Reference"])
            else:
                bad_writer.writerow([asv, best_sseqid, best_mismatch, "Other"])

if __name__ == "__main__":
    main()
