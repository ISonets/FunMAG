#!/usr/bin/env python3
"""
Calculate binning statistics: percentage of contigs binned by each method.
Usage: binned_mags_stats.py metabat_dir maxbin2_dir bin3c_dir contigs.fasta output.txt
"""

import sys
import os

def count_contigs_in_bins(bin_dir, extensions=None):
    """Count unique contig IDs across all bin files in a directory."""
    if extensions is None:
        extensions = ['.fa', '.fna', '.fasta']
    contigs = set()
    if not os.path.isdir(bin_dir):
        return 0
    for fname in os.listdir(bin_dir):
        if any(fname.endswith(ext) for ext in extensions):
            with open(os.path.join(bin_dir, fname)) as f:
                for line in f:
                    if line.startswith('>'):
                        contig_id = line[1:].strip().split()[0]
                        contigs.add(contig_id)
    return len(contigs)

def main():
    if len(sys.argv) != 6:
        print("Usage: binned_mags_stats.py <metabat_dir> <maxbin2_dir> <bin3c_dir> <contigs.fasta> <output.txt>",
              file=sys.stderr)
        sys.exit(1)

    metabat_dir = sys.argv[1]
    maxbin2_dir = sys.argv[2]
    bin3c_dir = sys.argv[3]
    contigs_file = sys.argv[4]
    out_file = sys.argv[5]

    # Count total contigs
    total_contigs = 0
    with open(contigs_file) as f:
        for line in f:
            if line.startswith('>'):
                total_contigs += 1

    if total_contigs == 0:
        with open(out_file, 'w') as out:
            out.write("Error: No contigs found in assembly.\n")
        sys.exit(1)

    metabat_n = count_contigs_in_bins(metabat_dir)
    maxbin2_n = count_contigs_in_bins(maxbin2_dir)
    bin3c_n = count_contigs_in_bins(bin3c_dir, extensions=['.fna'])

    with open(out_file, 'w') as out:
        out.write(f"Total contigs in assembly: {total_contigs}\n\n")
        out.write(f"MetaBAT2: {metabat_n:>6} contigs  ({100*metabat_n/total_contigs:5.1f}%)\n")
        out.write(f"MaxBin2:  {maxbin2_n:>6} contigs  ({100*maxbin2_n/total_contigs:5.1f}%)\n")
        out.write(f"bin3C:    {bin3c_n:>6} contigs  ({100*bin3c_n/total_contigs:5.1f}%)\n")

    print(f"Binning statistics written to {out_file}", file=sys.stderr)

if __name__ == '__main__':
    main()
