#!/bin/bash
# =============================================================================
# download_databases.sh – Download all reference databases
# Run AFTER install_envs.sh
# Adjust paths to match your config.yaml
# =============================================================================

set -e

# ---- Configuration (EDIT THESE) ----
DB_BASE="/path/to/databases"          # Root directory for all databases
THREADS=32                             # Threads for downloads and indexing

# Database paths (match config.yaml)
KRAKEN2_DB="${DB_BASE}/kraken2/standard"
BLAST_NR="${DB_BASE}/ncbi/nr"
KAIJU_DB="${DB_BASE}/kaiju/nr_euk"
KAIJU_NODES="${DB_BASE}/kaiju/nodes.dmp"
KRAKENUNIQ_DB="${DB_BASE}/krakenuniq/db"
BUSCO_DB="${DB_BASE}/busco"
HUMAN_DB="${DB_BASE}/kneaddata/human"
VIRSORTER2_DB="${DB_BASE}/virsorter2"

# ---- Color output ----
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  Database Downloads${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""

# ---- Create base directories ----
mkdir -p "${DB_BASE}"
mkdir -p "${BUSCO_DB}"
mkdir -p "${KRAKEN2_DB}"
mkdir -p "$(dirname ${BLAST_NR})"
mkdir -p "$(dirname ${KAIJU_DB})"
mkdir -p "${KRAKENUNIQ_DB}"
mkdir -p "${HUMAN_DB}"
mkdir -p "${VIRSORTER2_DB}"

# =============================================================================
# 1. BUSCO lineages (fungal)
# =============================================================================
echo -e "${GREEN}[1/7] Downloading BUSCO fungal lineages...${NC}"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate busco

export BUSCO_CONFIG_FILE="${BUSCO_DB}/config.ini"
mkdir -p "${BUSCO_DB}/lineages"

for lineage in fungi_odb10 ascomycota_odb10 basidiomycota_odb10 saccharomycetes_odb10; do
    if [ -d "${BUSCO_DB}/lineages/${lineage}" ]; then
        echo -e "${YELLOW}  ${lineage} already downloaded.${NC}"
    else
        echo "  Downloading ${lineage}..."
        busco --download ${lineage} --download_path "${BUSCO_DB}/lineages"
    fi
done

conda deactivate

# =============================================================================
# 2. Kraken2 standard database
# =============================================================================
echo ""
echo -e "${GREEN}[2/7] Building Kraken2 standard database...${NC}"
echo -e "${YELLOW}  WARNING: This is ~60 GB and takes several hours.${NC}"
echo -e "${YELLOW}  You can skip this if you already have a Kraken2 database.${NC}"
echo ""
read -p "  Download Kraken2 standard database? (y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    conda activate hic_mag
    if [ -f "${KRAKEN2_DB}/hash.k2d" ]; then
        echo -e "${YELLOW}  Kraken2 database exists.${NC}"
    else
        kraken2-build --standard --db "${KRAKEN2_DB}" --threads ${THREADS}
    fi
    conda deactivate
fi

# =============================================================================
# 3. BLAST nr database
# =============================================================================
echo ""
echo -e "${GREEN}[3/7] Downloading BLAST nr database...${NC}"
echo -e "${YELLOW}  WARNING: The nr database is ~250 GB compressed, ~350 GB uncompressed.${NC}"
echo -e "${YELLOW}  Total space needed: ~600 GB. This will take many hours.${NC}"
echo ""
read -p "  Download BLAST nr database? (y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    conda activate blast
    
    NR_DIR="$(dirname ${BLAST_NR})"
    mkdir -p "${NR_DIR}"
    
    if [ -f "${BLAST_NR}.ndb" ] || [ -f "${BLAST_NR}.00.nhd" ]; then
        echo -e "${YELLOW}  BLAST nr database already exists at ${BLAST_NR}.${NC}"
        echo -e "${YELLOW}  Delete it first to re-download.${NC}"
    else
        echo "  Step 1: Downloading nr FASTA files (this may take hours)..."
        # Download all nr FASTA files from NCBI
        wget --recursive --no-parent --accept 'nr.*.tar.gz' \
             --directory-prefix="${NR_DIR}/download" \
             https://ftp.ncbi.nlm.nih.gov/blast/db/ 2>&1 | tail -20
        
        echo "  Step 2: Extracting archives..."
        for tar_file in ${NR_DIR}/download/nr.*.tar.gz; do
            echo "    Extracting $(basename ${tar_file})..."
            tar -xzf "${tar_file}" -C "${NR_DIR}"
        done
        
        # Clean up downloaded archives to save space
        echo "  Step 3: Cleaning up downloaded archives..."
        rm -rf "${NR_DIR}/download"
        
        echo "  Step 4: Verifying database files..."
        if [ -f "${BLAST_NR}.ndb" ] || [ -f "${BLAST_NR}.00.nhd" ]; then
            echo -e "${GREEN}  BLAST nr database ready at ${BLAST_NR}${NC}"
        else
            echo -e "${RED}  ERROR: BLAST nr database files not found. Check download.${NC}"
            echo -e "${YELLOW}  Alternative: Use a pre-built nr database from:${NC}"
            echo -e "${YELLOW}  ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz${NC}"
            echo -e "${YELLOW}  Then run: makeblastdb -in nr -dbtype nucl -out nr${NC}"
        fi
    fi
    
    conda deactivate
fi

# =============================================================================
# 4. Kaiju nr_euk database
# =============================================================================
echo ""
echo -e "${GREEN}[4/7] Building Kaiju nr_euk database...${NC}"
echo -e "${YELLOW}  This is ~30 GB and takes several hours.${NC}"
echo ""
read -p "  Download Kaiju nr_euk database? (y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    conda activate kaiju
    if [ -f "${KAIJU_DB}.fmi" ]; then
        echo -e "${YELLOW}  Kaiju database exists.${NC}"
    else
        kaiju-makedb -s nr_euk -t ${THREADS}
        # Move to desired location if different from default
        mkdir -p "$(dirname ${KAIJU_DB})"
        if [ -f "kaiju_db_nr_euk.fmi" ]; then
            mv kaiju_db_nr_euk.fmi "${KAIJU_DB}.fmi"
            echo -e "${GREEN}  Kaiju database ready.${NC}"
        fi
    fi
    conda deactivate
fi

# =============================================================================
# 5. KrakenUniq database
# =============================================================================
echo ""
echo -e "${GREEN}[5/7] Building KrakenUniq database...${NC}"
echo -e "${YELLOW}  This requires an existing Kraken2/Bracken database as base.${NC}"
echo ""
read -p "  Build KrakenUniq database? (y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    conda activate krakenuniq
    if [ -f "${KRAKENUNIQ_DB}/database.kdb" ]; then
        echo -e "${YELLOW}  KrakenUniq database exists.${NC}"
    else
        # KrakenUniq can use standard Kraken2 db as source
        if [ -d "${KRAKEN2_DB}" ]; then
            krakenuniq-build --db "${KRAKENUNIQ_DB}" \
                             --kraken-db "${KRAKEN2_DB}" \
                             --threads ${THREADS}
        else
            echo -e "${RED}  Kraken2 database required first (step 2).${NC}"
        fi
    fi
    conda deactivate
fi

# =============================================================================
# 6. Human genome (for kneaddata decontamination)
# =============================================================================
echo ""
echo -e "${GREEN}[6/7] Downloading human reference genome (hg38)...${NC}"
echo -e "${YELLOW}  This is ~3 GB.${NC}"
echo ""
read -p "  Download human genome for decontamination? (y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    conda activate kneaddata
    if [ -d "${HUMAN_DB}" ] && [ "$(ls -A ${HUMAN_DB})" ]; then
        echo -e "${YELLOW}  Human database exists.${NC}"
    else
        kneaddata_database --download human_genome bowtie2 "${HUMAN_DB}"
        echo -e "${GREEN}  Human genome database ready.${NC}"
    fi
    conda deactivate
fi

# =============================================================================
# 7. VirSorter2 database
# =============================================================================
echo ""
echo -e "${GREEN}[7/7] Downloading VirSorter2 database...${NC}"
echo -e "${YELLOW}  This is ~2 GB.${NC}"
echo ""
read -p "  Download VirSorter2 database? (y/n): " -n 1 -r
echo
if [[ $REPLY =~ ^[Yy]$ ]]; then
    conda activate virsorter2
    if [ -f "${VIRSORTER2_DB}/checkv_database/hmms" ]; then
        echo -e "${YELLOW}  VirSorter2 database exists.${NC}"
    else
        virsorter setup --db-dir "${VIRSORTER2_DB}" -j ${THREADS}
    fi
    conda deactivate
fi

# =============================================================================
# Final summary
# =============================================================================
echo ""
echo -e "${GREEN}========================================${NC}"
echo -e "${GREEN}  Database download complete.${NC}"
echo -e "${GREEN}========================================${NC}"
echo ""
echo "Update config.yaml with these paths:"
echo "  kraken2_db:       \"${KRAKEN2_DB}\""
echo "  blast_nr_db:      \"${BLAST_NR}\""
echo "  kaiju_db:         \"${KAIJU_DB}\""
echo "  kaiju_nodes:      \"${KAIJU_NODES}\""
echo "  krakenuniq_db:    \"${KRAKENUNIQ_DB}\""
echo "  kneaddata_db:     \"${HUMAN_DB}\""
echo "  virsorter2_db:    \"${VIRSORTER2_DB}\""
echo ""
echo "BUSCO lineages are in: ${BUSCO_DB}/lineages"
echo "Set BUSCO_CONFIG_FILE to: ${BUSCO_DB}/config.ini"
