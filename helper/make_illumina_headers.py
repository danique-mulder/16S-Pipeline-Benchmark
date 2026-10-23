#!/usr/bin/env python3

# Rewrites FASTQ headers in Illumina format, so Natrix accepts the simulated reads in OTU mode.

import sys
import gzip

def write_fastq_with_new_headers(in_path, out_path, read_num, start_coord=10000):
    with gzip.open(in_path, "rt") as fin, gzip.open(out_path, "wt") as fout:
        idx = 0
        while True:
            header = fin.readline()
            if not header:
                break
            seq = fin.readline()
            plus = fin.readline()
            qual = fin.readline()
            if not (seq and plus and qual):
                break

            # Build unique cluster coordinates
            x = start_coord + idx
            y = start_coord + idx
            tile = 1101  # fixed tile number

            new_header = f"@INST:123:FC001:1:{tile}:{x}:{y} {read_num}:N:0:ACGTACGT\n"

            fout.write(new_header)
            fout.write(seq)
            fout.write(plus)
            fout.write(qual)

            idx += 1

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: make_illumina_headers.py <R1_in.fastq.gz> <R2_in.fastq.gz> <output_prefix>")
        sys.exit(1)

    r1_in = sys.argv[1]
    r2_in = sys.argv[2]
    out_prefix = sys.argv[3]

    write_fastq_with_new_headers(r1_in, f"{out_prefix}_R1.fastq.gz", read_num=1)
    write_fastq_with_new_headers(r2_in, f"{out_prefix}_R2.fastq.gz", read_num=2)
