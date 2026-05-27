#!/usr/bin/env python3
"""
Run BUSCO on a list of bins and filter by completeness/contamination.
Supports --auto-lineage for BUSCO 5.x.
Usage: python run_busco_filter.py --bin-list bins.lst --output results.tsv
       --filtered passing.lst --lineage auto --threads N
       --min-completeness 50 --max-contamination 10 [--mode comprehensive]
"""

import sys
import os
import argparse
import subprocess
import csv
import json
import tempfile

def parse_busco_json(json_file):
    """Extract completeness and contamination from BUSCO short_summary.json."""
    with open(json_file) as f:
        data = json.load(f)
    results = data.get('results', data)  # handle both BUSCO v4 and v5 formats
    single = results.get('single_copy', results.get('complete_single_copy', 0))
    duplicated = results.get('duplicated', results.get('complete_duplicated', 0))
    fragmented = results.get('fragmented', 0)
    missing = results.get('missing', 0)
    total = single + duplicated + fragmented + missing
    if total == 0:
        return 0.0, 0.0
    completeness = (single + duplicated) / total * 100
    contamination = duplicated / total * 100
    return round(completeness, 2), round(contamination, 2)

def run_busco_single(bin_path, output_dir, lineage, threads):
    """Run BUSCO on a single bin. Returns path to JSON file or None."""
    bin_name = os.path.basename(bin_path).replace('.fna', '').replace('.fa', '').replace('.fasta', '')
    busco_out = os.path.join(output_dir, bin_name)
    os.makedirs(busco_out, exist_ok=True)

    # Try requested lineage first
    if lineage == 'auto':
        cmd = (f"busco -i {bin_path} -o {bin_name} "
               f"--auto-lineage -m genome -c {threads} "
               f"--out_path {output_dir} --offline")
    else:
        cmd = (f"busco -i {bin_path} -o {bin_name} "
               f"-l {lineage} -m genome -c {threads} "
               f"--out_path {output_dir} --offline")

    try:
        subprocess.run(cmd, shell=True, check=True, capture_output=True, timeout=7200)
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired):
        # Fallback to fungi_odb10 if first attempt fails
        fallback_name = f"{bin_name}_fallback"
        fallback_out = os.path.join(output_dir, fallback_name)
        os.makedirs(fallback_out, exist_ok=True)
        cmd2 = (f"busco -i {bin_path} -o {fallback_name} "
                f"-l fungi_odb10 -m genome -c {threads} "
                f"--out_path {output_dir} --offline")
        try:
            subprocess.run(cmd2, shell=True, check=True, capture_output=True, timeout=7200)
            busco_out = fallback_out
        except:
            return None

    # Find the JSON file – check the specific output directory
    for root, dirs, files in os.walk(busco_out):
        for f in files:
            if f.startswith('short_summary') and f.endswith('.json'):
                return os.path.join(root, f)
    return None

def main():
    parser = argparse.ArgumentParser(description='Run BUSCO on bins and filter by quality.')
    parser.add_argument('--bin-list', required=True, help='File with paths to bin FASTA files')
    parser.add_argument('--output', required=True, help='Output TSV with BUSCO scores')
    parser.add_argument('--filtered', required=True, help='Output list of passing bins')
    parser.add_argument('--lineage', default='fungi_odb10', help='BUSCO lineage or "auto"')
    parser.add_argument('--threads', type=int, default=16, help='Threads per BUSCO run')
    parser.add_argument('--min-completeness', type=float, default=50.0)
    parser.add_argument('--max-contamination', type=float, default=10.0)
    parser.add_argument('--mode', default='comprehensive', choices=['comprehensive', 'dastool_only'])
    args = parser.parse_args()

    with open(args.bin_list) as f:
        bins = [line.strip() for line in f if line.strip()]

    if not bins:
        print("No bins to process.", file=sys.stderr)
        with open(args.output, 'w', newline='') as out:
            writer = csv.writer(out, delimiter='\t')
            writer.writerow(['bin', 'completeness', 'contamination', 'lineage', 'status'])
        with open(args.filtered, 'w') as f:
            pass
        return

    print(f"Running BUSCO on {len(bins)} bins...", file=sys.stderr)

    results = []
    passed = []
    busco_tmp = tempfile.mkdtemp(prefix='busco_tmp_')

    for i, bin_path in enumerate(bins):
        bin_name = os.path.basename(bin_path).replace('.fna', '').replace('.fa', '').replace('.fasta', '')
        print(f"  [{i+1}/{len(bins)}] {bin_name}...", file=sys.stderr, end=' ')

        json_file = run_busco_single(bin_path, busco_tmp, args.lineage, args.threads)

        if json_file and os.path.exists(json_file):
            comp, cont = parse_busco_json(json_file)
            status = 'PASS' if comp >= args.min_completeness and cont <= args.max_contamination else 'FAIL'
            lineage_used = os.path.basename(os.path.dirname(json_file))
            print(f"Completeness={comp:.1f}%, Contamination={cont:.1f}% -> {status}", file=sys.stderr)
        else:
            comp, cont = 0.0, 100.0
            status = 'FAIL'
            lineage_used = 'N/A'
            print(f"BUSCO failed -> FAIL", file=sys.stderr)

        results.append([bin_path, comp, cont, lineage_used, status])
        if status == 'PASS':
            passed.append(bin_path)

    # Write results
    with open(args.output, 'w', newline='') as out:
        writer = csv.writer(out, delimiter='\t')
        writer.writerow(['bin', 'completeness', 'contamination', 'lineage', 'status'])
        writer.writerows(results)

    # Write filtered list
    with open(args.filtered, 'w') as f:
        for p in passed:
            f.write(p + '\n')

    # Cleanup
    import shutil
    shutil.rmtree(busco_tmp, ignore_errors=True)

    print(f"\n{len(passed)}/{len(bins)} bins passed BUSCO filter", file=sys.stderr)

if __name__ == '__main__':
    main()
