#!/usr/bin/env python3
"""
Prepare contig-to-bin TSV files for DASTool from multiple binners.
Usage: python prep_tsv_for_dastool_multi.py
       metabat_lst.txt metabat metabat.tsv
       maxbin2_lst.txt maxbin2 maxbin2.tsv
       hic_lst.txt bin3c hic.tsv
"""

import sys
import os

def process_binner(lst_file, tag, out_tsv):
    """Convert a contig list file (one line per bin, contigs space-separated)
       into a contig-to-bin TSV file."""
    if not os.path.exists(lst_file):
        print(f"Warning: {lst_file} not found, creating empty {out_tsv}", file=sys.stderr)
        with open(out_tsv, 'w'):
            pass
        return

    with open(lst_file) as f, open(out_tsv, 'w') as out:
        for bin_id, line in enumerate(f, 1):
            contigs = line.strip().split()
            for contig in contigs:
                out.write(f"{contig}\t{tag}_{bin_id}\n")

def main():
    if len(sys.argv) < 4 or len(sys.argv) % 3 != 1:
        print("Usage: prep_tsv_for_dastool_multi.py <lst1> <tag1> <tsv1> [<lst2> <tag2> <tsv2> ...]",
              file=sys.stderr)
        sys.exit(1)

    args = sys.argv[1:]
    for i in range(0, len(args), 3):
        lst_file, tag, out_tsv = args[i], args[i+1], args[i+2]
        process_binner(lst_file, tag, out_tsv)

    print("DASTool TSV preparation complete.", file=sys.stderr)

if __name__ == '__main__':
    main()
