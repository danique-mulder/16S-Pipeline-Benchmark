
#!/usr/bin/env python3
import argparse
from Bio import SeqIO

def clean_genus_label(label):
    if label is None:
        return None
    label = label.strip()
    if label.startswith("g__"):
        label = label[3:]
    return label.strip()

def read_merged_taxonomy(file_path):
    """
    Read merged taxonomy file with header: SILVA_ID, exp, q2, idtaxa
    Clean genus labels.
    Returns dict: silvaID -> (exp, q2, idtaxa), all cleaned.
    """
    tax = {}
    with open(file_path) as f:
        header = f.readline()
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 4:
                continue
            silva_id = parts[0]
            exp = clean_genus_label(parts[1])
            q2 = clean_genus_label(parts[2])
            idtaxa = clean_genus_label(parts[3])
            tax[silva_id] = (exp, q2, idtaxa)
    return tax

def read_taxonomy(file_path):
    tax = {}
    with open(file_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                taxid = parts[0]
                tax_label = parts[1]
                tax[taxid] = tax_label
    return tax

def read_mapping(file_path):
    mapping = {}
    with open(file_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                silva_id, asvs = parts
                mapping[silva_id] = asvs
    return mapping

def write_filtered_count(file_path, input_count, valid_ids):
    with open(input_count) as f_in, open(file_path, "w") as f_out:
        header = f_in.readline().strip()
        f_out.write(header + "\n")
        for line in f_in:
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            sseqid = parts[0]
            if sseqid in valid_ids:
                f_out.write(line)

def write_filtered_taxonomy(file_path, taxonomy_dict, valid_ids):
    with open(file_path, 'w') as f:
        for sid in valid_ids:
            if sid in taxonomy_dict:
                f.write(f"{sid}\t{taxonomy_dict[sid]}\n")

def write_filtered_fasta(input_fasta, output_fasta, valid_ids):
    with open(output_fasta, 'w') as out_f:
        for record in SeqIO.parse(input_fasta, 'fasta'):
            if record.id in valid_ids:
                SeqIO.write(record, out_f, 'fasta')

def write_filtered_mapping(file_path, mapping, valid_ids):
    with open(file_path, 'w') as f:
        for sid in valid_ids:
            if sid in mapping:
                f.write(f"{sid}\t{mapping[sid]}\n")

def main(args):
    merged_tax = read_merged_taxonomy(args.merged_tax)
    taxonomy = read_taxonomy(args.taxonomy)
    mapping = read_mapping(args.mapping)

    merged_lines_out = []
    merged_lines_out.append("SILVA_ID\texp\tq2\tidtaxa\tall_match")

    valid_ids = set()
    for sid, (exp, q2, idtaxa) in merged_tax.items():
        if exp and q2 and idtaxa:
            all_match = (exp.lower() == q2.lower() == idtaxa.lower())
        else:
            all_match = False
        merged_lines_out.append(f"{sid}\t{exp}\t{q2}\t{idtaxa}\t{str(all_match)}")
        if all_match:
            valid_ids.add(sid)

    # write merged file with all_match column
    with open(args.out_merge, 'w') as f:
        f.write('\n'.join(merged_lines_out) + '\n')

    # filter outputs
    write_filtered_taxonomy(args.out_tax, taxonomy, valid_ids)
    write_filtered_count(args.out_count, args.count, valid_ids)
    write_filtered_fasta(args.fasta, args.out_fasta, valid_ids)
    write_filtered_mapping(args.out_mapping, mapping, valid_ids)

    print(f"Total SILVA IDs in merged file: {len(merged_tax)}")
    print(f"SILVA IDs with matching genus labels across all three: {len(valid_ids)}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter based on genus label agreement in merged taxonomy file.")
    parser.add_argument('--merged_tax', required=True, help="Merged taxonomy input file (with header: SILVA_ID exp q2 idtaxa)")
    parser.add_argument('--count', required=True, help="count_table_resolvable.tsv (sseqid qseqid taxonomy)")
    parser.add_argument('--taxonomy', required=True, help="Original taxonomy file (SILVA ID and taxonomy)")
    parser.add_argument('--fasta', required=True, help="FASTA file with sequences (IDs = SILVA IDs)")
    parser.add_argument('--mapping', required=True, help="Mapping file (silvaID<TAB>ASV1;ASV2;...)")
    parser.add_argument('--out_merge', required=True, help="Output merged taxonomy file with all_match column")
    parser.add_argument('--out_tax', required=True, help="Output filtered taxonomy file")
    parser.add_argument('--out_count', required=True, help="Output filtered count table (same format as count_table_resolvable.tsv)")
    parser.add_argument('--out_fasta', required=True, help="Output filtered FASTA file")
    parser.add_argument('--out_mapping', required=True, help="Output filtered mapping file")
    args = parser.parse_args()

    main(args)