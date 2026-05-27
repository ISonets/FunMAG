#!/usr/bin/env python3
"""
Run MAGPurify on fungal bins — safe modules only (tetra-freq, gc-content, coverage).
Skips phylo-markers and clade-markers (bacterial/archaeal databases).
Usage: python run_magpurify.py --bin-list bins.lst --output out_dir
       --bam alignment.bam --threads 16
"""

import sys
import os
import argparse
import subprocess

def main():
    parser = argparse.ArgumentParser(description='MAGPurify for fungal bins (safe modules).')
    parser.add_argument('--bin-list', required=True, help='File with paths to bins')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--bam', required=True, help='Sorted BAM file for coverage module')
    parser.add_argument('--threads', type=int, default=16)
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    with open(args.bin_list) as f:
        bins = [line.strip() for line in f if line.strip()]

    if not bins:
        print("No bins to process.", file=sys.stderr)
        return

    print(f"Running MAGPurify on {len(bins)} fungal bins...", file=sys.stderr)
    print("Using safe modules: tetra-freq, gc-content, coverage", file=sys.stderr)

    for i, bin_file in enumerate(bins):
        bin_name = os.path.basename(bin_file).replace('.fa', '').replace('.fna', '').replace('.fasta', '')
        print(f"  [{i+1}/{len(bins)}] {bin_name}", file=sys.stderr)

        out_bin = os.path.join(args.output, f"{bin_name}_clean.fa")
        success = True

        try:
            subprocess.run(f"magpurify tetra-freq {bin_file} {args.output}",
                           shell=True, capture_output=True, timeout=300)
        except Exception as e:
            print(f"    tetra-freq failed: {e}", file=sys.stderr)
            success = False

        try:
            subprocess.run(f"magpurify gc-content {bin_file} {args.output}",
                           shell=True, capture_output=True, timeout=300)
        except Exception as e:
            print(f"    gc-content failed: {e}", file=sys.stderr)

        if os.path.exists(args.bam) and os.path.getsize(args.bam) > 0:
            try:
                subprocess.run(f"magpurify coverage {bin_file} {args.output} --bam {args.bam}",
                               shell=True, capture_output=True, timeout=600)
            except Exception as e:
                print(f"    coverage failed: {e}", file=sys.stderr)

        if success:
            try:
                subprocess.run(f"magpurify clean-bin {bin_file} {out_bin} --threads {args.threads}",
                               shell=True, capture_output=True, timeout=600)
                if os.path.exists(out_bin):
                    print(f"    -> {out_bin}", file=sys.stderr)
                else:
                    raise RuntimeError("clean-bin produced no output")
            except Exception as e:
                print(f"    clean-bin failed: {e}", file=sys.stderr)
                success = False

        if not success or not os.path.exists(out_bin):
            print(f"    Copying original as fallback", file=sys.stderr)
            import shutil
            shutil.copy(bin_file, out_bin)

    print("MAGPurify complete.", file=sys.stderr)
    
if __name__ == '__main__':
    main()
