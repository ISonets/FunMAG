#!/usr/bin/env python3
"""
Generate two phased haplotype FASTA files from a phased VCF.
Uses bcftools consensus to apply phased variants to the reference.
"""

import sys
import os
import argparse
import subprocess
import tempfile

def main():
    parser = argparse.ArgumentParser(description='Generate phased haplotypes from VCF.')
    parser.add_argument('--vcf', required=True, help='Phased VCF file')
    parser.add_argument('--reference', required=True, help='Reference FASTA')
    parser.add_argument('--output1', required=True, help='Output haplotype 1 FASTA')
    parser.add_argument('--output2', required=True, help='Output haplotype 2 FASTA')
    args = parser.parse_args()

    # Check VCF is non-empty
    if not os.path.exists(args.vcf) or os.path.getsize(args.vcf) == 0:
        print("VCF empty, skipping haplotype generation", file=sys.stderr)
        # Copy reference as both haplotypes
        subprocess.run(f"cp {args.reference} {args.output1}", shell=True)
        subprocess.run(f"cp {args.reference} {args.output2}", shell=True)
        return

    print(f"Generating phased haplotypes from {args.vcf}...", file=sys.stderr)

    # Create haplotype 1 (0|1 alleles → take ref for 0, alt for 1)
    # For phased genotypes: 0|1 means haplotype 1 = ref, haplotype 2 = alt
    cmd1 = (f"bcftools consensus -H 1 -f {args.reference} {args.vcf} "
            f"2>/dev/null > {args.output1}")
    subprocess.run(cmd1, shell=True, check=True)

    # Create haplotype 2 (0|1 alleles → take alt for 1)
    cmd2 = (f"bcftools consensus -H 2 -f {args.reference} {args.vcf} "
            f"2>/dev/null > {args.output2}")
    subprocess.run(cmd2, shell=True, check=True)

    # Check if outputs were created
    for out in [args.output1, args.output2]:
        if os.path.exists(out) and os.path.getsize(out) > 0:
            print(f"  Created: {out} ({os.path.getsize(out)} bytes)", file=sys.stderr)
        else:
            print(f"  Warning: {out} is empty, using reference as fallback", file=sys.stderr)
            subprocess.run(f"cp {args.reference} {out}", shell=True)

    print("Phasing complete.", file=sys.stderr)


if __name__ == '__main__':
    main()
    
