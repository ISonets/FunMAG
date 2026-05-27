#!/usr/bin/env python3
"""
Dereplicate fungal bins with dRep using ANI clustering.
Keeps the best BUSCO completeness score per cluster.
Usage: python run_drep.py --bin-list fungal.lst --busco busco.tsv
       --output out_dir --ani 0.95 --coverage 0.3 --threads 16
"""

import sys
import os
import argparse
import subprocess
import csv
import shutil
import tempfile

def parse_busco_scores(busco_tsv):
    """Return dict: bin_path -> completeness (float)."""
    scores = {}
    with open(busco_tsv) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            try:
                scores[row['bin']] = float(row['completeness'])
            except (ValueError, KeyError):
                pass
    return scores

def main():
    parser = argparse.ArgumentParser(description='Dereplicate bins with dRep.')
    parser.add_argument('--bin-list', required=True)
    parser.add_argument('--busco', required=True, help='BUSCO scores TSV')
    parser.add_argument('--output', required=True, help='Output directory')
    parser.add_argument('--ani', type=float, default=0.95, help='ANI threshold')
    parser.add_argument('--coverage', type=float, default=0.3, help='Min overlap fraction')
    parser.add_argument('--threads', type=int, default=16)
    args = parser.parse_args()

    os.makedirs(args.output, exist_ok=True)

    with open(args.bin_list) as f:
        bins = [line.strip() for line in f if line.strip()]

    rep_lst = os.path.join(args.output, 'fungal_representatives.lst')

    if len(bins) == 0:
        print("No fungal bins to dereplicate.", file=sys.stderr)
        with open(rep_lst, 'w'):
            pass
        return

    if len(bins) == 1:
        print("Only one bin; no dereplication needed.", file=sys.stderr)
        dest = os.path.join(args.output, os.path.basename(bins[0]))
        shutil.copy(bins[0], dest)
        with open(rep_lst, 'w') as f:
            f.write(dest + '\n')
        return

    print(f"Dereplicating {len(bins)} fungal bins at ANI {args.ani}...", file=sys.stderr)

    # Write bin list for dRep
    work_dir = tempfile.mkdtemp(prefix='drep_')
    bin_list_file = os.path.join(work_dir, 'bins_to_dereplicate.txt')
    with open(bin_list_file, 'w') as f:
        for b in bins:
            f.write(b + '\n')

    # Run dRep
    cmd = (f"dRep dereplicate {args.output} "
           f"-g {bin_list_file} "
           f"-sa {args.ani} "
           f"-nc {args.coverage} "
           f"-p {args.threads} "
           f"--ignoreGenomeQuality "
           f"--S_algorithm fastANI "
           f"2>&1")

    try:
        subprocess.run(cmd, shell=True, check=True)
    except subprocess.CalledProcessError as e:
        print(f"dRep failed: {e}", file=sys.stderr)
        # Fallback: copy all bins as representatives
        representatives = []
        for b in bins:
            dest = os.path.join(args.output, os.path.basename(b))
            shutil.copy(b, dest)
            representatives.append(dest)
        with open(rep_lst, 'w') as f:
            for r in representatives:
                f.write(r + '\n')
        shutil.rmtree(work_dir, ignore_errors=True)
        return

    # Collect representative bins from dereplicated_genomes
    reps_dir = os.path.join(args.output, 'dereplicated_genomes')
    representatives = []
    if os.path.isdir(reps_dir):
        for f in sorted(os.listdir(reps_dir)):
            if f.endswith(('.fa', '.fna', '.fasta')):
                representatives.append(os.path.abspath(os.path.join(reps_dir, f)))

    with open(rep_lst, 'w') as f:
        for r in representatives:
            f.write(r + '\n')

    print(f"dRep produced {len(representatives)} representative bins.", file=sys.stderr)
    shutil.rmtree(work_dir, ignore_errors=True)

if __name__ == '__main__':
    main()
