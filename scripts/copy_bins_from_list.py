#!/usr/bin/env python3
"""
Copy bin files from a list to a destination directory.
Usage: python copy_bins_from_list.py bin_list.txt destination_dir/
"""

import sys
import os
import shutil

def main():
    if len(sys.argv) != 3:
        print("Usage: copy_bins_from_list.py <bin_list> <dest_dir>", file=sys.stderr)
        sys.exit(1)

    bin_list = sys.argv[1]
    dest_dir = sys.argv[2]

    if not os.path.exists(bin_list):
        print(f"Bin list file not found: {bin_list}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(dest_dir, exist_ok=True)

    copied = 0
    with open(bin_list) as f:
        for line in f:
            src = line.strip()
            if src and os.path.exists(src):
                shutil.copy(src, dest_dir)
                copied += 1
            elif src:
                print(f"Warning: file not found: {src}", file=sys.stderr)

    print(f"Copied {copied} bins to {dest_dir}", file=sys.stderr)

if __name__ == '__main__':
    main()
