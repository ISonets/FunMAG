#!/bin/bash
# =============================================================================
# Install all conda environments for the complete fungal MAG project.
# Covers: Core pipeline, Juicer, Hi‑C analysis, Post-review, Comparative genomics
# Run: bash install_envs.sh
# =============================================================================

set -e  # Exit on error
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENV_DIR="${SCRIPT_DIR}/envs"

# ---- Color output ----
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m' # No Color

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  Fungal MAG Project – Environment Setup${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""

# ---- Check conda ----
if ! command -v conda &> /dev/null; then
    echo -e "${RED}ERROR: conda not found. Install Miniconda/Anaconda first.${NC}"
    exit 1
fi

CONDA_VERSION=$(conda --version 2>&1)
echo -e "Conda version: ${CONDA_VERSION}"

# ---- Detect architecture ----
ARCH=$(uname -m)
OS=$(uname -s)
echo "System: ${OS} ${ARCH}"

# ---- Check for mamba (optional but recommended) ----
if command -v mamba &> /dev/null; then
    MAMBA_AVAILABLE=true
    echo -e "${GREEN}mamba detected – will use for faster installs.${NC}"
else
    MAMBA_AVAILABLE=false
    echo -e "${YELLOW}Tip: Install mamba for faster environment creation:${NC}"
    echo -e "${YELLOW}  conda install -n base -c conda-forge mamba${NC}"
fi
echo ""

# ---- Add channels globally ----
echo -e "${YELLOW}Configuring conda channels...${NC}"
conda config --add channels conda-forge 2>/dev/null || true
conda config --add channels bioconda 2>/dev/null || true
conda config --add channels defaults 2>/dev/null || true
conda config --set channel_priority flexible

# ---- Create environment function ----
create_env() {
    local env_name=$1
    local yaml_file="${ENV_DIR}/${env_name}.yaml"
    
    if [ ! -f "$yaml_file" ]; then
        echo -e "${RED}  ERROR: ${yaml_file} not found. Skipping.${NC}"
        return 1
    fi
    
    if conda env list | grep -q "^${env_name} "; then
        echo -e "${YELLOW}  Environment '${env_name}' already exists. Skipping.${NC}"
        return 0
    fi
    
    echo -e "${GREEN}  Creating '${env_name}'...${NC}"
    
    if $MAMBA_AVAILABLE; then
        mamba env create -f "$yaml_file" -q 2>&1 | tail -3
    else
        conda env create -f "$yaml_file" -q 2>&1 | tail -3
    fi
    
    echo -e "${GREEN}    Done.${NC}"
}

# =============================================================================
# Core Pipeline Environments
# =============================================================================
echo ""
echo -e "${CYAN}========================================${NC}"
echo -e "${CYAN}  Core MAG Pipeline Environments${NC}"
echo -e "${CYAN}========================================${NC}"

create_env "bioinf_v1"
create_env "hic_mag"
create_env "hicmag_py37"
create_env "bin3c"
create_env "maxbin2"
create_env "dastool"
create_env "busco"
create_env "blast"
create_env "kaiju"
create_env "drep"
create_env "metawrap"
create_env "magpurify"
create_env "virsorter2"
create_env "coverm"
create_env "qc3c"
create_env "kneaddata"
create_env "krakenuniq"
create_env "quast"
create_env "multiqc"

# =============================================================================
# HiCzin (pip install into its conda environment)
# =============================================================================
echo ""
echo -e "${CYAN}========================================${NC}"
echo -e "${CYAN}  Hi‑C Processing Environments${NC}"
echo -e "${CYAN}========================================${NC}"

if conda env list | grep -q "^hiczin "; then
    echo -e "${YELLOW}  HiCzin environment exists.${NC}"
else
    echo -e "${GREEN}  Creating HiCzin environment...${NC}"
    conda env create -f "${ENV_DIR}/hiczin.yaml" -q 2>&1 | tail -3
    CONDA_BASE=$(conda info --base)
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate hiczin
    pip install hiczin 2>&1 | tail -3 || echo -e "${RED}  HiCzin pip install failed – may need manual install${NC}"
    conda deactivate
    echo -e "${GREEN}    HiCzin installed.${NC}"
fi

# =============================================================================
# Juicer Environment
# =============================================================================
echo ""
echo -e "${CYAN}========================================${NC}"
echo -e "${CYAN}  Juicer Hi‑C Pipeline Environment${NC}"
echo -e "${CYAN}========================================${NC}"

create_env "juicer"

# Try to install Juicer + juicer_tools via conda
echo ""
echo -e "${YELLOW}Checking Juicer installation...${NC}"
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate juicer 2>/dev/null || true

if command -v juicer.sh &> /dev/null; then
    echo -e "${GREEN}  juicer.sh found in PATH.${NC}"
else
    echo -e "${YELLOW}  juicer.sh not in conda. Installing from GitHub...${NC}"
    bash "${SCRIPT_DIR}/scripts/install_juicer.sh" "${HOME}/tools/juicer"
fi

if command -v juicer_tools.jar &> /dev/null; then
    echo -e "${GREEN}  juicer_tools.jar found in PATH.${NC}"
else
    echo -e "${YELLOW}  juicer_tools.jar not found. Will search at runtime.${NC}"
fi

conda deactivate 2>/dev/null || true

# =============================================================================
# Comparative Genomics Environment (anvi'o)
# =============================================================================
echo ""
echo -e "${CYAN}========================================${NC}"
echo -e "${CYAN}  Comparative Genomics (anvi'o) Environment${NC}"
echo -e "${CYAN}========================================${NC}"

create_env "anvio"

# Setup KEGG KOfams for anvi'o
echo ""
echo -e "${YELLOW}Setting up KEGG KOfams for anvi'o...${NC}"
echo -e "${YELLOW}  This requires kegg-db.tar.gz in the project root.${NC}"
echo -e "${YELLOW}  Download from: https://www.genome.jp/kegg-bin/get_htext?ko00001${NC}"
echo ""

if [ -f "${SCRIPT_DIR}/kegg-db.tar.gz" ]; then
    source "$(conda info --base)/etc/profile.d/conda.sh"
    conda activate anvio
    anvi-setup-kegg-kofams --kegg-archive "${SCRIPT_DIR}/kegg-db.tar.gz" 2>&1 | tail -5
    conda deactivate
    echo -e "${GREEN}  KEGG KOfams setup complete.${NC}"
else
    echo -e "${YELLOW}  kegg-db.tar.gz not found. Skipping KEGG setup.${NC}"
    echo -e "${YELLOW}  Run manually after downloading:${NC}"
    echo -e "${YELLOW}    conda activate anvio${NC}"
    echo -e "${YELLOW}    anvi-setup-kegg-kofams --kegg-archive kegg-db.tar.gz${NC}"
fi

# =============================================================================
# Verify installations
# =============================================================================
echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  Verifying installations...${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""

CONDA_BASE=$(conda info --base)
source "${CONDA_BASE}/etc/profile.d/conda.sh"

verify_tool() {
    local env=$1
    local tool=$2
    local version_flag=${3:---version}
    
    conda activate "$env" 2>/dev/null
    if command -v "$tool" &> /dev/null; then
        version=$($tool $version_flag 2>&1 | head -1 | cut -c1-80)
        echo -e "  ${GREEN}✓${NC} ${env}/${tool}: ${version}"
    else
        echo -e "  ${RED}✗${NC} ${env}/${tool}: NOT FOUND"
    fi
    conda deactivate 2>/dev/null
}

echo -e "${CYAN}Core pipeline:${NC}"
verify_tool "bioinf_v1"     "fastqc"      "--version"
verify_tool "bioinf_v1"     "samtools"    "--version"
verify_tool "hic_mag"       "bwa"
verify_tool "hic_mag"       "metaphlan"   "--version"
verify_tool "hic_mag"       "kraken2"     "--version"
verify_tool "hicmag_py37"   "spades.py"   "--version"
verify_tool "bin3c"         "bin3C.py"    "--version"
verify_tool "busco"         "busco"       "--version"
verify_tool "blast"         "blastn"      "-version"
verify_tool "kaiju"         "kaiju"       "-h" 2>/dev/null | head -1
verify_tool "drep"          "dRep"        "--version"
verify_tool "coverm"        "coverm"      "--version"
verify_tool "kneaddata"     "kneaddata"   "--version"
verify_tool "krakenuniq"    "krakenuniq"  "--version"
verify_tool "quast"         "quast.py"    "--version"

echo ""
echo -e "${CYAN}Hi‑C pipeline:${NC}"
verify_tool "juicer"        "bwa"
verify_tool "juicer"        "samtools"    "--version"
# juicer.sh and juicer_tools.jar checked separately above

echo ""
echo -e "${CYAN}Comparative genomics:${NC}"
verify_tool "anvio"         "anvi-pan-genome" "--version" 2>/dev/null || echo -e "  ${YELLOW}⚠${NC} anvi'o – check manually"
verify_tool "anvio"         "augustus"    "--version"
verify_tool "anvio"         "muscle"      "-version"
verify_tool "anvio"         "iqtree"      "--version"

echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  Installation complete.${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo -e "${CYAN}Next steps:${NC}"
echo ""
echo -e "  ${GREEN}1.${NC} Update database paths in config.yaml"
echo -e "  ${GREEN}2.${NC} Download required databases:"
echo -e "      bash download_databases.sh"
echo ""
echo -e "  ${GREEN}3.${NC} Run the core pipeline:"
echo -e "      snakemake --cores 60 --use-conda"
echo ""
echo -e "  ${GREEN}4.${NC} Run Juicer Hi‑C pipeline:"
echo -e "      snakemake --cores 32 --use-conda -s Snakefile_juicer"
echo ""
echo -e "  ${GREEN}5.${NC} Run Hi‑C analysis:"
echo -e "      snakemake --cores 16 --use-conda -s Snakefile_hic_analysis"
echo ""
echo -e "  ${GREEN}6.${NC} After Juicebox review:"
echo -e "      bash post_review.sh"
echo ""
echo -e "  ${GREEN}7.${NC} Run comparative genomics:"
echo -e "      Place FASTA files in fastas_comp_genomics/"
echo -e "      snakemake --cores 32 --use-conda -s Snakefile_comp_genomics"
echo ""
echo -e "${YELLOW}Note: Some tools (MaxBin2, DAS_Tool, metaWRAP) use Perl/Ruby scripts${NC}"
echo -e "${YELLOW}and may need their environment activated to verify properly.${NC}"
echo ""
echo -e "${YELLOW}KEGG database for anvi'o:${NC}"
echo -e "${YELLOW}  Download from: https://www.genome.jp/kegg-bin/get_htext?ko00001${NC}"
echo -e "${YELLOW}  Place as: kegg-db.tar.gz in project root${NC}"
