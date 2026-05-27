#!/usr/bin/env python3
"""
Merge Arrowhead domain calls from multiple resolutions into a single BEDPE file.
"""
import sys
import os
import glob
import argparse

def merge_arrowhead_output(input_dir, output_bedpe):
    """Find all Arrowhead output files and merge into one BEDPE."""
    all_domains = []
    
    # Arrowhead outputs files named like: 5000_blocks.bedpe, 10000_blocks.bedpe, etc.
    pattern = os.path.join(input_dir, "*_blocks.bedpe")
    for bedpe_file in glob.glob(pattern):
        if os.path.exists(bedpe_file) and os.path.getsize(bedpe_file) > 0:
            with open(bedpe_file) as f:
                for line in f:
                    if line.strip() and not line.startswith('#'):
                        all_domains.append(line.strip())
    
    with open(output_bedpe, 'w') as out:
        out.write("#chr1\tstart1\tend1\tchr2\tstart2\tend2\tcolor\tscore\n")
        for domain in all_domains:
            out.write(domain + '\n')
    
    print(f"Merged {len(all_domains)} domains from {input_dir}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('input_dir', help='Directory with Arrowhead output')
    parser.add_argument('output_bedpe', help='Output merged BEDPE file')
    args = parser.parse_args()
    merge_arrowhead_output(args.input_dir, args.output_bedpe)

if __name__ == '__main__':
    main()
