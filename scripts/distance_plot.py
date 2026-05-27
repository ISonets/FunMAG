#!/usr/bin/env python3
"""
Generate distance-dependent contact frequency plot from .hic file.
Requires juicer_tools and matplotlib.
"""
import sys
import os
import subprocess
import tempfile
import argparse
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

def extract_contact_vs_distance(hic_file, juicer_tools, resolution=10000):
    """Extract contact probability vs distance using juicer_tools."""
    tmp_dir = tempfile.mkdtemp()
    try:
        cmd = (f"java -Xmx8g -jar {juicer_tools} dump observed KR "
               f"{hic_file} chr1 chr1 BP {resolution} "
               f"{tmp_dir}/contacts.txt")
        subprocess.run(cmd, shell=True, check=True, capture_output=True)
        
        distances = []
        contacts = []
        
        with open(f"{tmp_dir}/contacts.txt") as f:
            for line in f:
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        pos1, pos2, count = int(parts[0]), int(parts[1]), float(parts[2])
                        dist = abs(pos2 - pos1)
                        if dist > 0:
                            distances.append(dist)
                            contacts.append(count)
        
        if not distances:
            return np.array([]), np.array([])
        
        # Bin by distance
        max_dist = max(distances)
        bins = np.logspace(3, np.log10(max_dist), 50)
        bin_indices = np.digitize(distances, bins)
        
        mean_contacts = []
        for i in range(1, len(bins)):
            mask = bin_indices == i
            if mask.any():
                mean_contacts.append(np.mean(np.array(contacts)[mask]))
            else:
                mean_contacts.append(0)
        
        bin_centers = (bins[:-1] + bins[1:]) / 2
        return bin_centers, np.array(mean_contacts)
    finally:
        import shutil
        shutil.rmtree(tmp_dir, ignore_errors=True)

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--hic', required=True)
    parser.add_argument('--output', required=True)
    parser.add_argument('--juicer-tools', required=True)
    parser.add_argument('--sample', default='unknown')
    args = parser.parse_args()
    
    distances, contacts = extract_contact_vs_distance(args.hic, args.juicer_tools)
    
    if len(distances) == 0:
        print("No contacts found", file=sys.stderr)
        plt.figure(figsize=(8, 6))
        plt.text(0.5, 0.5, 'No data', ha='center', va='center')
        plt.savefig(args.output, dpi=150)
        return
    
    plt.figure(figsize=(10, 6))
    plt.loglog(distances, contacts, 'b-', linewidth=2, alpha=0.7)
    plt.xlabel('Genomic Distance (bp)', fontsize=12)
    plt.ylabel('Contact Probability', fontsize=12)
    plt.title(f'Distance vs Contact Frequency – {args.sample}', fontsize=14)
    plt.grid(True, alpha=0.3)
    plt.tight_layout()
    plt.savefig(args.output, dpi=150)
    print(f"Plot saved to {args.output}", file=sys.stderr)

if __name__ == '__main__':
    main()
