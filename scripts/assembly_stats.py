#!/usr/bin/env python3
"""
Calculate assembly statistics: N50, L50, N90, total length, GC content, etc.
"""
import sys
import argparse
from Bio import SeqIO
import numpy as np

def calculate_assembly_stats(fasta_file):
    """Calculate comprehensive assembly statistics."""
    lengths = []
    gc_count = 0
    total_bases = 0
    num_contigs = 0
    num_ns = 0

    for record in SeqIO.parse(fasta_file, 'fasta'):
        seq = str(record.seq).upper()
        contig_len = len(record.seq)
        lengths.append(contig_len)
        total_bases += contig_len
        num_contigs += 1
        gc_count += seq.count('G') + seq.count('C')
        num_ns += seq.count('N')

    if not lengths:
        return None

    lengths = sorted(lengths, reverse=True)
    half_total = total_bases / 2

    # N50, L50
    cumsum = 0
    n50 = 0
    l50 = 0
    for i, l in enumerate(lengths):
        cumsum += l
        if cumsum >= half_total:
            n50 = l
            l50 = i + 1
            break

    # N90, L90
    ninety_percent = total_bases * 0.9
    cumsum = 0
    n90 = 0
    l90 = 0
    for i, l in enumerate(lengths):
        cumsum += l
        if cumsum >= ninety_percent:
            n90 = l
            l90 = i + 1
            break

    stats = {
        'total_length': total_bases,
        'num_contigs': num_contigs,
        'max_contig': lengths[0],
        'min_contig': lengths[-1],
        'mean_contig': np.mean(lengths),
        'median_contig': np.median(lengths),
        'n50': n50,
        'l50': l50,
        'n90': n90,
        'l90': l90,
        'gc_content': (gc_count / total_bases * 100) if total_bases > 0 else 0,
        'num_ns': num_ns,
        'percent_ns': (num_ns / total_bases * 100) if total_bases > 0 else 0
    }
    return stats

def main():
    parser = argparse.ArgumentParser(description='Calculate assembly statistics.')
    parser.add_argument('--fasta', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--sample', default='unknown')
    args = parser.parse_args()

    stats = calculate_assembly_stats(args.fasta)

    if stats is None:
        with open(args.output, 'w') as f:
            f.write("No contigs found.\n")
        return

    with open(args.output, 'w') as f:
        f.write(f"Assembly Statistics for {args.sample}\n")
        f.write(f"{'='*50}\n\n")
        f.write(f"Total length:      {stats['total_length']:>15,} bp\n")
        f.write(f"Number of contigs: {stats['num_contigs']:>15,}\n")
        f.write(f"Longest contig:    {stats['max_contig']:>15,} bp\n")
        f.write(f"Shortest contig:   {stats['min_contig']:>15,} bp\n")
        f.write(f"Mean contig size:  {stats['mean_contig']:>15,.1f} bp\n")
        f.write(f"Median contig:     {stats['median_contig']:>15,.1f} bp\n")
        f.write(f"N50:               {stats['n50']:>15,} bp\n")
        f.write(f"L50:               {stats['l50']:>15,}\n")
        f.write(f"N90:               {stats['n90']:>15,} bp\n")
        f.write(f"L90:               {stats['l90']:>15,}\n")
        f.write(f"GC content:        {stats['gc_content']:>15.2f}%\n")
        f.write(f"N's in assembly:   {stats['num_ns']:>15,} ({stats['percent_ns']:.2f}%)\n")

    print(f"Assembly stats written to {args.output}", file=sys.stderr)

if __name__ == '__main__':
    main()
