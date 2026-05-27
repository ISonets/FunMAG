#!/usr/bin/env python3
"""
Calculate library complexity statistics from Juicer merged_nodups.txt.
Computes PBC1, PBC2, NRF (Non-Redundant Fraction), and unique read pairs.
Based on ENCODE Hi‑C library complexity metrics.
"""

import sys
import argparse
from collections import defaultdict

def parse_merged_nodups(filepath):
    """Parse Juicer merged_nodups.txt and count read pairs."""
    total_pairs = 0
    unique_positions = set()  # (chr1, pos1, strand1, chr2, pos2, strand2)
    unique_reads = set()      # read_name

    with open(filepath) as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split()
            if len(parts) < 7:
                continue

            # Format: read_name str1 chr1 pos1 frag1 str2 chr2 pos2 frag2 mapq1 mapq2 ...
            read_name = parts[0]
            str1 = parts[1]
            chr1 = parts[2]
            pos1 = int(parts[3])
            str2 = parts[6]
            chr2 = parts[7]
            pos2 = int(parts[8])

            total_pairs += 1
            unique_reads.add(read_name)
            unique_positions.add((chr1, pos1, str1, chr2, pos2, str2))

    return total_pairs, len(unique_reads), len(unique_positions)


def calculate_metrics(total_pairs, distinct_reads, distinct_positions):
    """Calculate PBC1, PBC2, NRF."""
    # PBC1 = distinct_reads / total_pairs
    pbc1 = distinct_reads / total_pairs if total_pairs > 0 else 0

    # PBC2 = distinct_positions / total_pairs
    pbc2 = distinct_positions / total_pairs if total_pairs > 0 else 0

    # NRF (Non-Redundant Fraction) = distinct_reads / distinct_positions
    nrf = distinct_reads / distinct_positions if distinct_positions > 0 else 0

    return pbc1, pbc2, nrf


def main():
    parser = argparse.ArgumentParser(description='Calculate Hi‑C library complexity.')
    parser.add_argument('--merged-nodups', required=True, help='Path to merged_nodups.txt')
    parser.add_argument('--output', required=True, help='Output metrics file')
    parser.add_argument('--sample', default='unknown', help='Sample name')
    args = parser.parse_args()

    print(f"Processing {args.merged_nodups}...", file=sys.stderr)

    total, distinct_reads, distinct_positions = parse_merged_nodups(args.merged_nodups)

    if total == 0:
        print("Warning: No read pairs found in merged_nodups.txt", file=sys.stderr)
        with open(args.output, 'w') as out:
            out.write("sample\ttotal_pairs\tdistinct_reads\tdistinct_positions\tPBC1\tPBC2\tNRF\n")
            out.write(f"{args.sample}\t0\t0\t0\t0\t0\t0\n")
        return

    pbc1, pbc2, nrf = calculate_metrics(total, distinct_reads, distinct_positions)

    with open(args.output, 'w') as out:
        out.write("sample\ttotal_pairs\tdistinct_reads\tdistinct_positions\tPBC1\tPBC2\tNRF\n")
        out.write(f"{args.sample}\t{total}\t{distinct_reads}\t{distinct_positions}\t"
                  f"{pbc1:.6f}\t{pbc2:.6f}\t{nrf:.6f}\n")

    print(f"  Total pairs: {total}", file=sys.stderr)
    print(f"  Distinct reads: {distinct_reads}", file=sys.stderr)
    print(f"  Distinct positions: {distinct_positions}", file=sys.stderr)
    print(f"  PBC1: {pbc1:.4f}", file=sys.stderr)
    print(f"  PBC2: {pbc2:.4f}", file=sys.stderr)
    print(f"  NRF: {nrf:.4f}", file=sys.stderr)
    print(f"Metrics written to {args.output}", file=sys.stderr)


if __name__ == '__main__':
    main()
