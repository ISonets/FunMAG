#!/usr/bin/env python3
"""
Compute median coverage across samples for MaxBin2 bins.
Usage: python median_cov_maxbin2.py depth.tsv maxbin2_bins_dir output.tsv
"""

import sys
import os
import re
import pandas as pd
import numpy as np

def main():
    if len(sys.argv) != 4:
        print("Usage: median_cov_maxbin2.py <depth.tsv> <bins_dir> <output.tsv>", file=sys.stderr)
        sys.exit(1)

    depth_file = sys.argv[1]
    bins_dir = sys.argv[2]
    out_file = sys.argv[3]

    depth = pd.read_csv(depth_file, sep='\t')
    cols_to_drop = [c for c in depth.columns if c.endswith('var')]
    depth = depth.drop(columns=cols_to_drop)
    if depth.shape[1] > 3:
        depth = depth.drop(depth.columns[[1, 2]], axis=1)

    bin_files = [f for f in os.listdir(bins_dir) if f.endswith(('.fa', '.fasta'))]
    if not bin_files:
        with open(out_file, 'w') as f:
            f.write("bin\tmedian_coverage\n")
        return

    pattern = re.compile(r'^>(\S+)')
    results = {}

    for fname in sorted(bin_files):
        bin_path = os.path.join(bins_dir, fname)
        contigs = set()
        with open(bin_path) as f:
            for line in f:
                if line.startswith('>'):
                    match = pattern.match(line)
                    if match:
                        contigs.add(match.group(1))

        if not contigs:
            continue

        sub = depth[depth['contigName'].isin(contigs)]
        if sub.empty or sub.shape[1] < 2:
            continue

        medians = sub.iloc[:, 1:].median(axis=0, numeric_only=True)
        bin_name = fname.replace('.fa', '').replace('.fasta', '')
        results[bin_name] = medians

    if not results:
        with open(out_file, 'w') as f:
            f.write("bin\tmedian_coverage\n")
        return

    out_df = pd.DataFrame(results).T
    out_df.index.name = 'bin'
    out_df.to_csv(out_file, sep='\t')

if __name__ == '__main__':
    main()
