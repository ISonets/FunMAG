#!/usr/bin/env python3
"""
Generate a summary HTML table from Kaiju output for all bins in a sample.
Complements the Krona chart with a sortable table.
"""

import sys
import os
import argparse
from collections import Counter, defaultdict

def parse_kaiju_output(kaiju_file):
    """
    Parse Kaiju output: bin_path<TAB>genus1;genus2;...
    Returns dict: bin_name -> list of genera
    """
    bins_taxonomy = {}
    with open(kaiju_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            bin_path = parts[0]
            bin_name = os.path.basename(bin_path).replace('.fa', '').replace('.fna', '').replace('.fasta', '')
            genera = [g.strip() for g in parts[1].split(';') if g.strip()]
            bins_taxonomy[bin_name] = genera
    return bins_taxonomy

def generate_html_table(bins_taxonomy, output_file):
    """Generate an HTML table with taxonomy summary."""
    # Count top genera across all bins
    genus_counter = Counter()
    bin_count = len(bins_taxonomy)
    
    for genera in bins_taxonomy.values():
        for g in genera:
            genus_counter[g] += 1
    
    top_genera = genus_counter.most_common(20)
    
    html = """
    <!DOCTYPE html>
    <html>
    <head>
        <title>Kaiju Taxonomy Summary</title>
        <style>
            body { font-family: Arial, sans-serif; margin: 20px; }
            table { border-collapse: collapse; width: 100%; }
            th, td { border: 1px solid #ddd; padding: 8px; text-align: left; }
            th { background-color: #4CAF50; color: white; }
            tr:nth-child(even) { background-color: #f2f2f2; }
            .summary { margin: 20px 0; padding: 15px; background: #e7f3ff; border-radius: 5px; }
        </style>
    </head>
    <body>
        <h1>Kaiju Taxonomy Summary</h1>
        <div class="summary">
            <p><strong>Total bins analyzed:</strong> BIN_COUNT</p>
            <p><strong>Unique genera detected:</strong> GENUS_COUNT</p>
        </div>
        <h2>Top 20 Most Abundant Genera</h2>
        <table>
            <tr><th>Genus</th><th>Number of Bins</th><th>Percentage</th></tr>
            TABLE_ROWS
        </table>
        <h2>Per-Bin Taxonomy</h2>
        <table>
            <tr><th>Bin Name</th><th>Assigned Genera</th></tr>
            BIN_ROWS
        </table>
    </body>
    </html>
    """
    
    table_rows = ""
    for genus, count in top_genera:
        pct = (count / bin_count * 100) if bin_count > 0 else 0
        table_rows += f"<tr><td>{genus}</td><td>{count}</td><td>{pct:.1f}%</td></tr>\n"
    
    bin_rows = ""
    for bin_name, genera in sorted(bins_taxonomy.items()):
        genera_str = ", ".join(genera[:10])  # Top 10 per bin
        if len(genera) > 10:
            genera_str += f" ... (+{len(genera)-10} more)"
        bin_rows += f"<tr><td>{bin_name}</td><td>{genera_str}</td></tr>\n"
    
    html = html.replace("BIN_COUNT", str(bin_count))
    html = html.replace("GENUS_COUNT", str(len(genus_counter)))
    html = html.replace("TABLE_ROWS", table_rows)
    html = html.replace("BIN_ROWS", bin_rows)
    
    with open(output_file, 'w') as f:
        f.write(html)

def main():
    parser = argparse.ArgumentParser(description='Generate Kaiju HTML summary table.')
    parser.add_argument('--kaiju', required=True, help='Kaiju output file')
    parser.add_argument('--output', required=True, help='Output HTML file')
    args = parser.parse_args()
    
    if not os.path.exists(args.kaiju):
        with open(args.output, 'w') as f:
            f.write("<html><body><p>No Kaiju data available.</p></body></html>")
        return
    
    bins_taxonomy = parse_kaiju_output(args.kaiju)
    generate_html_table(bins_taxonomy, args.output)
    print(f"HTML report written to {args.output}", file=sys.stderr)

if __name__ == '__main__':
    main()
