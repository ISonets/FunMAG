# 🍄 Fungal MAG Reconstruction & Comparative Genomics Pipeline

[![Snakemake](https://img.shields.io/badge/Snakemake-≥7.0-brightgreen)](https://snakemake.github.io/)
[![Python](https://img.shields.io/badge/Python-3.7|3.9-blue)](https://www.python.org/)
[![Conda](https://img.shields.io/badge/Conda-environments-green)](https://docs.conda.io/)
[![License](https://img.shields.io/badge/License-MIT-yellow)](LICENSE)

A comprehensive, modular Snakemake-based pipeline for **metagenome-assembled genome (MAG) reconstruction**, **Hi‑C contact map generation and analysis**, and **comparative genomics** with a focus on **fungal genomes from various environments**.

---

## 📋 Table of Contents

- [Overview](#-overview)
- [Pipeline Architecture](#-pipeline-architecture)
- [Installation](#-installation)
- [Configuration](#-configuration)
- [Usage](#-usage)
- [Output Structure](#-output-structure)
- [Resource Requirements](#-resource-requirements)
- [Pipeline Flow](#-pipeline-flow)
- [Methods Description](#-methods-description)
- [Troubleshooting](#-troubleshooting)
- [Citation](#-citation)
- [License](#-license)

---

## 🔬 Overview

This pipeline processes paired-end **WGS (Whole Genome Shotgun)** and **Hi‑C (chromosome conformation capture)** sequencing data to:

1. **Reconstruct fungal MAGs** from metagenomic assemblies using multiple binning approaches
2. **Generate Hi‑C contact maps** for each fungal MAG using Juicer
3. **Analyze Hi‑C maps** — normalization, compartments, domains, loops, APA
4. **Manually review and correct** assemblies in Juicebox GUI
5. **Perform comparative genomics** — pangenome analysis, phylogenomics, ANI, functional enrichment, and metabolism estimation using anvi'o

The pipeline is designed for **beverage fermentation metagenomes** (wine, beer, sake, kombucha) but is applicable to any fungal-dominated metagenomic dataset.

---

## 🏗 Pipeline Architecture

The project consists of **5 interconnected Snakefiles**:

| # | Snakefile | Purpose |
|---|-----------|---------|
| 1 | `snake_core_fungal_bins` | **Core MAG reconstruction** — assembly → binning → BUSCO → taxonomy → dRep → MetaWRAP → MAGPurify |
| 2 | `juicer` | **Juicer Hi‑C mapping** — contact map generation for each fungal MAG |
| 3 | `hic_map_analysis` | **Hi‑C map analysis** — normalization, compartments, domains, loops, APA |
| 4 | `juicebox_post_review` | **Post-Juicebox review** — corrected FASTA, new .hic, QUAST, comparison |
| 5 | `comp_genomics` | **Comparative genomics** — anvi'o pangenome, phylogenomics, ANI, enrichment |

All pipelines share configuration files and conda environments, ensuring reproducibility.

---

## 💾 Installation

### Prerequisites

- **Conda** (Miniconda or Anaconda)
- **Snakemake** ≥ 7.0
- **Git**
- **Java** ≥ 11 (for Juicer tools)
- **Perl** ≥ 5.32 (for Juicer statistics)

### Quick Install

```bash
# Clone the repository
git clone https://github.com/your-username/fungal-mag-pipeline.git
cd fungal-mag-pipeline

# Install all conda environments (~30 minutes)
bash install_envs.sh

# Download required databases (interactive)
bash download_databases.sh

### Database Requirements

| Database | Size | Used By |
|----------|------|---------|
| BUSCO lineages (fungi_odb10, ascomycota_odb10, etc.) | ~2 GB | Core pipeline |
| Kraken2 standard database | ~60 GB | Core pipeline |
| BLAST nr database | ~350 GB | Core pipeline |
| Kaiju nr_euk database | ~30 GB | Core pipeline |
| KrakenUniq microbial database | ~20 GB | Core pipeline |
| Human genome (hg38, for kneaddata) | ~3 GB | Core pipeline |
| VirSorter2 database | ~2 GB | Core pipeline |
| KEGG KOfams (for anvi'o) | ~1 GB | Comparative genomics |

See `download_databases.sh` for automated download instructions.

---

## ⚙ Configuration

### Core Pipeline (`config.yaml`)

```yaml
samples: [""] #e.g. B32
data_dir: "/path/to/raw/data"
n_threads: 32
spades_mem: 30  # Adjust for your RAM
restriction_enzyme: "" #e.g. HpaII
busco_lineage: "auto"
drep_ani_threshold: 0.95  # Species level
```

### Juicer Pipeline (`config_juicer.yaml`)

```yaml
samples: [""] #e.g. B32
restriction_enzyme: "" #e.g. HpaII
hic_resolutions: "5000,10000,25000,50000,100000,250000,500000,1000000"
run_3d_dna: true
run_phasing: false  # Enable only for diploid fungi
```

### Comparative Genomics (`config_comp_genomics.yaml`)

```yaml
project_name: "" #e.g. Brettanomyces_bruxellensis
genome_dir: "fastas_comp_genomics"
augustus_species: "" #e.g. "pichia_stipitis"
mcl_inflation: 10
```

**Important:** Update all database paths in the config files before running.

---

## 🚀 Usage

### 1. Prepare Input Data

```bash
# Organize raw FASTQ files
mkdir -p data/B32
ln -s /path/to/B32_WGS_R1.fastq.gz data/B32/
ln -s /path/to/B32_WGS_R2.fastq.gz data/B32/
ln -s /path/to/B32_HiC_R1.fastq.gz data/B32/
ln -s /path/to/B32_HiC_R2.fastq.gz data/B32/
```

### 2. Run Core Pipeline

```bash
# Dry-run first to verify
snakemake --cores 32 --use-conda -n -s snake_core_fungal_bins

# Full execution
snakemake --cores 32 --use-conda 2>&1 -s snake_core_fungal_bins | tee core_pipeline.log
```

### 3. Run Juicer Hi‑C Mapping

```bash
snakemake --cores 32 --use-conda -s juicer 2>&1 | tee juicer.log
```

### 4. Run Hi‑C Map Analysis

```bash
snakemake --cores 32 --use-conda -s hic_map_analysis 2>&1 | tee hic_analysis.log
```

### 5. Manual Review in Juicebox

1. Open `juicer/{sample}/{mag_name}/aligned/inter_30.hic` in Juicebox GUI
2. Import assembly from `juicer/{sample}/{mag_name}/3d-dna/`
3. Manually correct misassemblies, inversions, translocations
4. Export corrected assembly as `genome.review.assembly`

### 6. Run Post-Review Pipeline

```bash
snakemake --cores 32 --use-conda -s juicer_post_review 2>&1 | tee juicer_post_review.log
```

### 7. Prepare for Comparative Genomics

```bash
# Copy corrected assemblies + reference genomes
mkdir -p fastas_comp_genomics
cp post_review/B32/mag_name/genome_corrected.fasta fastas_comp_genomics/Brettanomyces_B32.fasta
# Download reference genomes from NCBI
wget -O fastas_comp_genomics/CBS_2499.fasta https://...
```

### 8. Run Comparative Genomics

```bash
snakemake --cores 32 --use-conda -s comp_genomics 2>&1 | tee comp_genomics.log
```

---

## 📂 Output Structure

```
results/{sample}/
├── fungal_bins_final.lst              # Final fungal MAGs
├── fungal_bins_magpurify/             # Cleaned MAGs (FASTA)
├── quast/report.html                  # Assembly quality report
├── coverm_fungal_abundance.tsv        # MAG abundance across samples
├── virsorter2/                        # Viral sequence predictions
│
juicer/{sample}/{mag_name}/
├── aligned/
│   ├── merged_nodups.txt              # Deduplicated Hi‑C contacts
│   └── inter_30.hic                   # Hi‑C contact map
├── 3d-dna/contigs.fasta              # 3D-DNA scaffolded assembly
│
hic_analysis/{sample}/{mag_name}/
├── normalized/                        # KR, VC, VC_SQRT normalized maps
├── compartments/eigenvectors.txt      # A/B compartments
├── domains/arrowhead_domains.bedpe    # Contact domains
├── loops/hiccups_loops.bedpe          # Chromatin loops
├── apa/apa_results.txt                # Aggregate Peak Analysis
└── qc/distance_plot.png              # QC: distance vs contacts
│
post_review/{sample}/{mag_name}/
├── genome_corrected.fasta             # Manually corrected assembly
├── hic/inter_30_corrected.hic        # New .hic from corrected FASTA
└── comparison_report.html             # Original vs corrected comparison
│
anvio/{PROJECT}/
├── {PROJECT}-PAN.db                   # Pangenome database
├── phylogenomic-tree.txt              # IQ-TREE phylogenomic tree
├── ANI/ANI_percentage_identity.txt    # PyANI similarity matrix
├── functional_enrichment.txt          # Functional enrichment (GLM)
└── metabolism_estimation.txt          # KEGG module completion
```

---

## 📊 Resource Requirements

### Tested Hardware

| Component | Minimum | Recommended |
|-----------|---------|-------------|
| CPU | 16 threads | 32+ threads |
| RAM | 32 GB | 64 GB (for SPAdes) |
| Storage | 500 GB SSD | 1+ TB NVMe SSD |
| OS | Linux/macOS | Ubuntu 20.04+ |

### Performance Estimates (per sample, 57M WGS + 46M Hi‑C reads)

| Pipeline | Time (32 threads, 32 GB RAM) | Peak RAM | Disk Usage |
|----------|------------------------------|----------|------------|
| Core | 40-80 hours | 30 GB (SPAdes) | 150-250 GB |
| Juicer | 4-8 hours per MAG | 16 GB | 6-18 GB per MAG |
| Hi‑C Analysis | 2-4 hours per MAG | 16 GB | 0.5-2.5 GB per MAG |
| Post-Review | 1-3 hours per MAG | 16 GB | 0.2-1 GB per MAG |
| Comparative Genomics | 12-24 hours total | 32 GB | 10-40 GB |

---

## 🔄 Pipeline Flow

```
Raw Reads (WGS + Hi‑C)
    │
    ├── Decontamination (kneaddata, hg38)
    │
    ├── Quality Control (FastQC, SeqPrep, Trimmomatic, qc3C, bbduk)
    │
    ├── Assembly (metaSPAdes / MEGAHIT)
    │
    ├── Binning (MetaBAT2, MaxBin2, bin3C) → DASTool consensus
    │
    ├── BUSCO Quality Assessment → Taxonomy (BLAST nr, Kaiju, MetaPhlAn, MiCoP)
    │
    ├── Fungal Bin Selection (consensus) → dRep Dereplication
    │
    ├── MetaWRAP Reassembly → MAGPurify Cleaning → Final Fungal MAGs
    │
    ├── Juicer Hi‑C Mapping → .hic generation → 3D-DNA Scaffolding
    │
    ├── Hi‑C Analysis (KR/VC normalization, compartments, domains, loops, APA)
    │
    ├── Juicebox Manual Review → Post-Review Corrections
    │
    └── Comparative Genomics (anvi'o pangenome, phylogenomics, ANI, enrichment)
```


## 🔧 Troubleshooting

### Common Issues

| Problem | Solution |
|---------|----------|
| **SPAdes OOM (Out of Memory)** | Set `spades_mem: 30` or switch to `use_megahit: true` |
| **BLAST nr database too large** | Use `-max_target_seqs 10` (already in pipeline) or switch to `blastn -task megablast` |
| **BUSCO fails with `--auto-lineage`** | Falls back to `fungi_odb10` automatically; check offline data |
| **Juicer `-S early` fails** | Ensure restriction enzyme matches library preparation; check `generate_site_positions.py` output |
| **Anvi'o KEGG setup fails** | Download `kegg-db.tar.gz` manually from KEGG and place in project root |
| **Conda environment conflicts** | Use `mamba` instead of `conda` for faster resolution; recreate environments with `--conda-cleanup-envs` |
| **Disk space exhausted** | Clean intermediate files: `snakemake --delete-all-output` then rerun with `--cores` |
| **MAGPurify fails on fungal bins** | Pipeline uses only safe modules (`tetra-freq`, `gc-content`, `coverage`); skip `phylo-markers` |

### Logs

All rules generate logs and benchmarks in `logs/`. Check:
- `logs/{rule}/{sample}.log` — stdout/stderr
- `logs/{rule}/{sample}.benchmark.txt` — runtime, RAM usage

---

## 📄 Citation

If you use this pipeline in your research, please cite:

```
[Sonets IV]. (2026). Fungal MAG Reconstruction & Comparative Genomics Pipeline.
GitHub: https://github.com/your-username/FunMAG
```

Additionally, please cite the key tools used:

- **SPAdes:** Prjibelski et al. (2020) *Current Protocols in Bioinformatics*
- **MetaBAT2:** Kang et al. (2019) *PeerJ*
- **MaxBin2:** Wu et al. (2016) *Bioinformatics*
- **bin3C:** DeMaere & Darling (2019) *mSystems*
- **DAS Tool:** Sieber et al. (2018) *Nature Microbiology*
- **BUSCO:** Manni et al. (2021) *Molecular Biology and Evolution*
- **dRep:** Olm et al. (2017) *ISME Journal*
- **MetaWRAP:** Uritskiy et al. (2018) *Microbiome*
- **MAGPurify:** Nayfach et al. (2019) *Nature Biotechnology*
- **Juicer:** Durand et al. (2016) *Cell Systems*
- **3D-DNA:** Dudchenko et al. (2017) *Science*
- **anvi'o:** Eren et al. (2021) *Nature Microbiology*
- **IQ-TREE:** Minh et al. (2020) *Molecular Biology and Evolution*

---

## 🤝 Contributing

Contributions are welcome! Please:

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/amazing-feature`)
3. Commit your changes (`git commit -m 'Add amazing feature'`)
4. Push to the branch (`git push origin feature/amazing-feature`)
5. Open a Pull Request

Please report bugs and feature requests via [GitHub Issues](https://github.com/your-username/fungal-mag-pipeline/issues).

---

## 📜 License

This project is licensed under the MIT License — see the [LICENSE](LICENSE) file for details.

---

## 🙏 Acknowledgments

- **ENCODE Consortium** for Hi‑C processing standards
- **Meren Lab** for anvi'o and its extensive documentation
- **Aiden Lab** for Juicer, 3D-DNA, and Juicebox
- All open-source tool developers whose software made this pipeline possible

---

## 📧 Contact

For questions, suggestions, or collaboration:

- **GitHub:** (https://github.com/ISonets)

---

<p align="center">
    <b>🍄 Happy fungal genome hunting! 🧬</b>
</p>
