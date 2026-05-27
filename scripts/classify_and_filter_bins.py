#!/usr/bin/env python3
"""
Select fungal bins by consensus of taxonomy tools (BLAST nr, Kaiju, MetaPhlAn, MiCoP).
Intersects with BUSCO-passed bins. UNITE removed — nr BLAST is sufficient.
Inputs passed via Snakemake 'script' directive.
"""

import os
import csv
import shutil
from collections import Counter

# ---------- PARSERS ----------

def load_blast(csv_file):
    """Identify fungal bins from BLAST nr results by keyword matching.
    Format: bin_path<TAB>blast_hit_1<TAB>blast_hit_2... (paste -s output)
    """
    fungal = set()
    fungal_keywords = [
        'fungi', 'ascomycota', 'basidiomycota', 'mucoromycota',
        'chytridiomycota', 'saccharomycetales', 'eurotiales',
        'hypocreales', 'glomus', 'rhizophagus', 'saccharomyces',
        'aspergillus', 'penicillium', 'candida', 'fusarium',
        'kluyveromyces', 'brettanomyces', 'zygosaccharomyces',
        'torulaspora', 'hanseniaspora', 'metschnikowia',
        'pichia', 'ogataea', 'yarrowia', 'alternaria', 'botrytis',
        'cladosporium', 'colletotrichum', 'magnaporthe'
    ]
    with open(csv_file) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            # First field (before first tab) is the bin path
            fields = line.split('\t')
            bin_path = fields[0]
            # Check all remaining fields for fungal keywords
            desc = ' '.join(fields[1:]).lower()
            if any(kw in desc for kw in fungal_keywords):
                fungal.add(bin_path)
    return fungal

def load_kaiju(kaiju_file):
    """Parse Kaiju output; check if genera are known fungal genera."""
    fungal_genera = {
        'saccharomyces', 'aspergillus', 'candida', 'penicillium',
        'fusarium', 'trichoderma', 'neurospora', 'schizosaccharomyces',
        'cryptococcus', 'ustilago', 'puccinia', 'agaricus',
        'coprinopsis', 'laccaria', 'trametes', 'pleurotus',
        'kluyveromyces', 'brettanomyces', 'dekkera', 'zygosaccharomyces',
        'lachancea', 'torulaspora', 'hanseniaspora', 'metschnikowia',
        'pichia', 'ogataea', 'yarrowia', 'alternaria', 'botrytis',
        'cladosporium', 'colletotrichum', 'magnaporthe', 'podospora'
    }
    fungal = set()
    with open(kaiju_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) < 2:
                continue
            bin_path = parts[0]
            genera_str = parts[1]
            for g in genera_str.split(';'):
                if g.strip().lower() in fungal_genera:
                    fungal.add(bin_path)
                    break
    return fungal


def load_metaphlan(csv_file):
    """Parse MetaPhlAn output; check for fungal kingdom/phylum assignment."""
    fungal = set()
    with open(csv_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            tax = row.get('taxonomy', '')
            if any(x in tax for x in ['k__Fungi', 'p__Ascomycota', 'p__Basidiomycota',
                                       'p__Mucoromycota', 'p__Chytridiomycota',
                                       'p__Zoopagomycota', 'p__Blastocladiomycota',
                                       's__Saccharomyces']):
                fungal.add(row['bin'])
    return fungal


def load_micop(csv_file, min_frac=0.5):
    """Parse MiCoP output; bin is fungal if fungal_reads/total_reads > threshold."""
    fungal = set()
    with open(csv_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                frac = float(row.get('fraction', 0))
                if frac >= min_frac:
                    fungal.add(row['bin'])
            except (ValueError, KeyError):
                pass
    return fungal


def load_krakenuniq(report_file):
    """
    Parse KrakenUniq report for fungal bins.
    Placeholder — implement based on actual KrakenUniq output format.
    """
    return set()


# ---------- MAIN ----------

def main():
    snakemake = globals().get('snakemake', None)
    if snakemake is None:
        import sys
        print("This script must be run via Snakemake's 'script' directive.", file=sys.stderr)
        sys.exit(1)

    blast_csv   = snakemake.input.get('blast', '')
    kaiju_out   = snakemake.input.get('kaiju', '')
    mpa_csv     = snakemake.input.get('metaphlan', '')
    micop_csv   = snakemake.input.get('micop', '')
    busco_file  = snakemake.input.get('busco_filtered', '')
    kraken_report = snakemake.input.get('dastool_kraken', '')

    min_frac = snakemake.params.min_frac
    out_dir  = snakemake.params.out_dir

    # Load BUSCO-passed bins
    busco_pass = set()
    if busco_file and os.path.exists(busco_file):
        with open(busco_file) as f:
            busco_pass = {line.strip() for line in f if line.strip()}

    if not busco_pass:
        print("No bins passed BUSCO filtering. Creating empty output.", file=sys.stderr)
        os.makedirs(out_dir, exist_ok=True)
        with open(snakemake.output.lst, 'w'):
            pass
        return

    # Collect fungal sets from each classifier
    fungal_sets = []
    classifiers_used = 0

    if blast_csv and os.path.exists(blast_csv):
        fungal_sets.append(load_blast(blast_csv))
        classifiers_used += 1
    if kaiju_out and os.path.exists(kaiju_out):
        fungal_sets.append(load_kaiju(kaiju_out))
        classifiers_used += 1
    if mpa_csv and os.path.exists(mpa_csv):
        fungal_sets.append(load_metaphlan(mpa_csv))
        classifiers_used += 1
    if micop_csv and os.path.exists(micop_csv):
        fungal_sets.append(load_micop(micop_csv, min_frac))
        classifiers_used += 1
    if kraken_report and os.path.exists(kraken_report):
        fungal_sets.append(load_krakenuniq(kraken_report))
        classifiers_used += 1

    if classifiers_used == 0:
        print("No taxonomy tools were run. Creating empty output.", file=sys.stderr)
        os.makedirs(out_dir, exist_ok=True)
        with open(snakemake.output.lst, 'w'):
            pass
        return

    # Majority vote: bin is fungal if >= half of the tools agree
    all_bins = set().union(*fungal_sets) if fungal_sets else set()
    counter = Counter()
    for s in fungal_sets:
        for b in s:
            counter[b] += 1

    threshold = max(1, classifiers_used // 2)
    fungal_bins_consensus = {b for b, cnt in counter.items() if cnt >= threshold}

    # Intersect with BUSCO-passed bins
    final_fungal = fungal_bins_consensus & busco_pass

    print(f"Classifiers used: {classifiers_used}, threshold: {threshold}", file=sys.stderr)
    print(f"Bins flagged fungal by consensus: {len(fungal_bins_consensus)}", file=sys.stderr)
    print(f"After BUSCO filter: {len(final_fungal)}", file=sys.stderr)

    # Copy to output directory
    os.makedirs(out_dir, exist_ok=True)
    final_bins = []
    for b in sorted(final_fungal):
        if os.path.exists(b):
            dest = os.path.join(out_dir, os.path.basename(b))
            shutil.copy(b, dest)
            final_bins.append(dest)

    with open(snakemake.output.lst, 'w') as f:
        for b in final_bins:
            f.write(b + '\n')

    print(f"Copied {len(final_bins)} fungal bins to {out_dir}", file=sys.stderr)


if __name__ == '__main__':
    main()
