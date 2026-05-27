#!/usr/bin/env python3
"""
Filter FASTA file to contigs >= min_length bp.
Usage: python filter_contigs_by_length.py input.fasta output.fasta --min-length N
"""

import sys
from Bio import SeqIO

def main():
    if len(sys.argv) != 5 or sys.argv[3] != '--min-length':
        print("Usage: filter_contigs_by_length.py input.fasta output.fasta --min-length N", file=sys.stderr)
        sys.exit(1)

    in_fasta = sys.argv[1]
    out_fasta = sys.argv[2]
    min_len = int(sys.argv[4])

    kept = 0
    total = 0
    with open(out_fasta, 'w') as out:
        for rec in SeqIO.parse(in_fasta, 'fasta'):
            total += 1
            if len(rec) >= min_len:
                SeqIO.write(rec, out, 'fasta')
                kept += 1

    print(f"Kept {kept}/{total} contigs >= {min_len} bp", file=sys.stderr)

if __name__ == '__main__':
    main()
