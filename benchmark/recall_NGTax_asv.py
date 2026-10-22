#!/usr/bin/env python3
"""
pair_NGTax_asv.py

Inputs:
  blast_forward.tsv  - BLAST outfmt 6 for forward reads
  blast_reverse.tsv  - BLAST outfmt 6 for reverse reads

For each query (qseqid), the script returns:

  1. ASVs with 100% match in forward AND reverse.
     If >0 found → use these ; separated

  2. Otherwise, ASVs where (mismatches_forward + mismatches_reverse) == 1.

Output: TSV table with:
  query   ASV_0_mismatch	ASV_1_mismatch
"""

import sys
import csv
from collections import defaultdict

def parse_blast(path):
    """Parses BLAST outfmt 6. Returns dict: hits[qid][sid] = mismatch_count."""
    hits = defaultdict(lambda: defaultdict(list))
    with open(path) as fh:
        rdr = csv.reader(fh, delimiter="\t")
        for row in rdr:
            (
                qid, sid, pident, mismatch, length
            ) = row

            mismatches = int(mismatch)

            # mismatches = qlen - matches = qlen - (pident% of qlen)

            hits[qid][sid].append(mismatches)

    # Convert list of mismatches per ASV into their minimum mismatch (best match)
    best = defaultdict(dict)
    for qid in hits:
        for sid in hits[qid]:
            best[qid][sid] = min(hits[qid][sid])

    return best


def main():
    if len(sys.argv) != 3:
        print("Usage: python pair_NGTax_asv.py forward.tsv reverse.tsv")
        sys.exit(1)

    fwd_file = sys.argv[1]
    rev_file = sys.argv[2]

    fwd = parse_blast(fwd_file)
    rev = parse_blast(rev_file)

    all_queries = set(fwd.keys()) | set(rev.keys())

    print("sequence.ID\tASV_0_mismatch\tASV_1_mismatch")

    for q in sorted(all_queries):
        strict = []
        fallback = []

        asvs = set(fwd.get(q, {}).keys()) | set(rev.get(q, {}).keys())

        for a in asvs:
            f_mis = fwd.get(q, {}).get(a, 999)   # if missing treat as bad match
            r_mis = rev.get(q, {}).get(a, 999)
            total = f_mis + r_mis

            # Rule 1: both exact
            if f_mis == 0 and r_mis == 0:
                strict.append(a)
            # Rule 2: total mismatches ≤ 1
            elif total <= 1:
                fallback.append(a)

        # Priority: strict > fallback
        if strict:
            print(f"{q}\t{','.join(strict)}\tNA")
        else:
            if fallback:
                print(f"{q}\tNA\t{','.join(fallback)}")
            else:
                print(f"{q}\tNA\tNA")


if __name__ == "__main__":
    main()
