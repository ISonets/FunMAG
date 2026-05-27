# Snakefile_comp_genomics – anvi'o Pangenome & Comparative Genomics
# Final version with header cleaning fix and PyANI integration.
# Takes finalized fungal genomes and performs pangenome analysis,
# phylogenomics, ANI, functional enrichment, and metabolism estimation.

import os
import glob

configfile: "config_comp_genomics.yaml"

# ======================================================================
# Input: all .fasta files in fastas_comp_genomics/
# ======================================================================

GENOME_DIR = config.get("genome_dir", "fastas_comp_genomics")
GENOMES = glob.glob(f"{GENOME_DIR}/*.fasta")
GENOME_NAMES = [os.path.basename(g).replace('.fasta', '') for g in GENOMES]

# Auto-detect number of genomes if not set in config
N_GENOMES = config.get("n_genomes", 0)
if N_GENOMES == 0:
    N_GENOMES = len(GENOME_NAMES)

# Project name for pangenome
PROJECT_NAME = config.get("project_name", "fungal_pangenome")

# ======================================================================
# Final targets
# ======================================================================

rule all:
    input:
        f"anvio/{PROJECT_NAME}/{PROJECT_NAME}-PAN.db",
        f"anvio/{PROJECT_NAME}_GENOMES.db",
        f"anvio/{PROJECT_NAME}/phylogenomic-tree.txt",
        f"anvio/{PROJECT_NAME}/ANI/ANI_percentage_identity.txt",
        f"anvio/{PROJECT_NAME}/functional_enrichment.txt",
        f"anvio/{PROJECT_NAME}/metabolism_estimation.txt",
        f"anvio/{PROJECT_NAME}/pangenome_summary.done",

# ======================================================================
# 1. Reformat FASTA files (anvi'o requirements)
# ======================================================================

rule reformat_fasta:
    """
    Reformat FASTA for anvi'o: remove non-ATGC characters,
    simplify names, filter contigs < 1000 bp.
    """
    input: f"{GENOME_DIR}/{{genome}}.fasta"
    output:
        fixed = "anvio/fixed_fastas/{genome}_fixed.fasta",
        done = touch("anvio/fixed_fastas/{genome}_reformat.done")
    benchmark: "logs/anvio/reformat/{genome}.benchmark.txt"
    log: "logs/anvio/reformat/{genome}.log"
    conda: "envs/anvio.yaml"
    shell:
        "anvi-script-reformat-fasta {input} "
        "-o {output.fixed} "
        "-l 1000 "
        "--simplify-names "
        "--seq-type NT "
        "2>&1 | tee {log}"

# ======================================================================
# 2. Gene prediction with AUGUSTUS
# ======================================================================

rule augustus_gene_prediction:
    """
    Predict genes using AUGUSTUS with species-specific model.
    Output GFF3 format.
    """
    input:
        fixed = "anvio/fixed_fastas/{genome}_fixed.fasta",
        reformat_done = "anvio/fixed_fastas/{genome}_reformat.done"
    output:
        gff = "anvio/augustus/{genome}_ann.gff",
        done = touch("anvio/augustus/{genome}_augustus.done")
    benchmark: "logs/anvio/augustus/{genome}.benchmark.txt"
    log: "logs/anvio/augustus/{genome}.log"
    threads: config.get("augustus_threads", 8)
    params:
        species = config.get("augustus_species", "pichia_stipitis")
    conda: "envs/anvio.yaml"
    shell:
        "augustus --species={params.species} "
        "--gff3=on "
        "--softmasking=1 "
        "--genemodel=complete "
        "{input.fixed} "
        "> {output.gff} 2>&1 | tee {log}"

# ======================================================================
# 3. Clean GFF (remove header lines that break anvi'o)
# ======================================================================

rule clean_gff:
    """
    Remove comment lines from GFF that break anvi'o parsing.
    """
    input:
        gff = "anvio/augustus/{genome}_ann.gff",
        augustus_done = "anvio/augustus/{genome}_augustus.done"
    output:
        cleaned = "anvio/augustus/{genome}_ann_cleaned.gff",
        done = touch("anvio/augustus/{genome}_clean.done")
    benchmark: "logs/anvio/clean_gff/{genome}.benchmark.txt"
    log: "logs/anvio/clean_gff/{genome}.log"
    shell:
        "sed '/^#/d' {input.gff} > {output.cleaned} 2>&1 | tee {log}"

# ======================================================================
# 4. Convert GFF to anvi'o external gene calls
# ======================================================================

rule gff_to_external_calls:
    """
    Convert AUGUSTUS GFF to anvi'o external gene calls format.
    """
    input:
        gff = "anvio/augustus/{genome}_ann_cleaned.gff",
        clean_done = "anvio/augustus/{genome}_clean.done"
    output:
        calls = "anvio/augustus/{genome}_ann.txt",
        done = touch("anvio/augustus/{genome}_calls.done")
    benchmark: "logs/anvio/external_calls/{genome}.benchmark.txt"
    log: "logs/anvio/external_calls/{genome}.log"
    conda: "envs/anvio.yaml"
    shell:
        "anvi-script-augustus-output-to-external-gene-calls "
        "-i {input.gff} "
        "-o {output.calls} "
        "2>&1 | tee {log}"

# ======================================================================
# 5. Create contigs database
# ======================================================================

rule create_contigs_db:
    """
    Create anvi'o contigs database with external gene calls.
    """
    input:
        fixed = "anvio/fixed_fastas/{genome}_fixed.fasta",
        calls = "anvio/augustus/{genome}_ann.txt",
        calls_done = "anvio/augustus/{genome}_calls.done"
    output:
        db = "anvio/databases/{genome}.db",
        done = touch("anvio/databases/{genome}_db.done")
    benchmark: "logs/anvio/contigs_db/{genome}.benchmark.txt"
    log: "logs/anvio/contigs_db/{genome}.log"
    conda: "envs/anvio.yaml"
    shell:
        "anvi-gen-contigs-database "
        "-f {input.fixed} "
        "-o {output.db} "
        "--external-gene-calls {input.calls} "
        "--project-name \"{config[project_name]}\" "
        "2>&1 | tee {log}"

# ======================================================================
# 6. Run HMMs (single-copy core genes)
# ======================================================================

rule run_hmms:
    """
    Run HMM profiles on contigs database to identify single-copy genes.
    """
    input:
        db = "anvio/databases/{genome}.db",
        db_done = "anvio/databases/{genome}_db.done"
    output:
        done = touch("anvio/databases/{genome}_hmms.done")
    benchmark: "logs/anvio/hmms/{genome}.benchmark.txt"
    log: "logs/anvio/hmms/{genome}.log"
    threads: config.get("hmm_threads", 8)
    conda: "envs/anvio.yaml"
    shell:
        "anvi-run-hmms "
        "-c {input.db} "
        "--num-threads {threads} "
        "2>&1 | tee {log}"

# ======================================================================
# 7. Setup KEGG KOfams (run once)
# ======================================================================

rule setup_kegg:
    """
    Download and setup KEGG KOfams database for functional annotation.
    Runs once for the project.
    """
    output:
        done = touch("anvio/kegg_setup.done")
    benchmark: "logs/anvio/kegg_setup.benchmark.txt"
    log: "logs/anvio/kegg_setup.log"
    conda: "envs/anvio.yaml"
    shell:
        "anvi-setup-kegg-kofams "
        "--kegg-archive kegg-db.tar.gz "
        "2>&1 | tee {log} && "
        "touch {output.done}"

# ======================================================================
# 8. Run KEGG KOfams annotation
# ======================================================================

rule run_kegg:
    """
    Annotate genes with KEGG KOfams.
    """
    input:
        db = "anvio/databases/{genome}.db",
        hmms_done = "anvio/databases/{genome}_hmms.done",
        kegg_setup = "anvio/kegg_setup.done"
    output:
        done = touch("anvio/databases/{genome}_kegg.done")
    benchmark: "logs/anvio/kegg/{genome}.benchmark.txt"
    log: "logs/anvio/kegg/{genome}.log"
    threads: config.get("kegg_threads", 8)
    conda: "envs/anvio.yaml"
    shell:
        "anvi-run-kegg-kofams "
        "-c {input.db} "
        "--num-threads {threads} "
        "2>&1 | tee {log}"

# ======================================================================
# 9. Run Pfam annotation
# ======================================================================

rule run_pfam:
    """
    Annotate genes with Pfam domains using HMMER.
    """
    input:
        db = "anvio/databases/{genome}.db",
        hmms_done = "anvio/databases/{genome}_hmms.done"
    output:
        done = touch("anvio/databases/{genome}_pfam.done")
    benchmark: "logs/anvio/pfam/{genome}.benchmark.txt"
    log: "logs/anvio/pfam/{genome}.log"
    threads: config.get("pfam_threads", 8)
    conda: "envs/anvio.yaml"
    shell:
        "anvi-run-pfams "
        "-c {input.db} "
        "--num-threads {threads} "
        "2>&1 | tee {log}"

# ======================================================================
# 10. Create external-genomes.txt
# ======================================================================

rule create_external_genomes:
    """
    Generate external-genomes.txt file listing all genomes.
    """
    input:
        dbs = expand("anvio/databases/{genome}.db", genome=GENOME_NAMES),
        kegg = expand("anvio/databases/{genome}_kegg.done", genome=GENOME_NAMES)
    output:
        external = "anvio/external-genomes.txt",
        done = touch("anvio/external_genomes.done")
    benchmark: "logs/anvio/external_genomes.benchmark.txt"
    log: "logs/anvio/external_genomes.log"
    run:
        with open(output.external, 'w') as out:
            out.write("name\tcontigs_db_path\n")
            for genome in GENOME_NAMES:
                db_path = os.path.abspath(f"anvio/databases/{genome}.db")
                out.write(f"{genome}\t{db_path}\n")

# ======================================================================
# 11. Create genome storage
# ======================================================================

rule create_genome_storage:
    """
    Create anvi'o genome storage from external-genomes.txt.
    """
    input:
        external = "anvio/external-genomes.txt",
        external_done = "anvio/external_genomes.done"
    output:
        storage = f"anvio/{PROJECT_NAME}_GENOMES.db",
        done = touch(f"anvio/{PROJECT_NAME}_genome_storage.done")
    benchmark: "logs/anvio/genome_storage.benchmark.txt"
    log: "logs/anvio/genome_storage.log"
    conda: "envs/anvio.yaml"
    params:
        gene_caller = config.get("gene_caller", "AUGUSTUS")
    shell:
        "anvi-gen-genomes-storage "
        "-e {input.external} "
        "-o {output.storage} "
        "--gene-caller {params.gene_caller} "
        "2>&1 | tee {log}"

# ======================================================================
# 12. Pangenome analysis
# ======================================================================

rule pangenome_analysis:
    """
    Run pangenome analysis with MCL inflation = 10.
    Uses DIAMOND for sequence alignment in sensitive mode.
    """
    input:
        storage = f"anvio/{PROJECT_NAME}_GENOMES.db",
        storage_done = f"anvio/{PROJECT_NAME}_genome_storage.done"
    output:
        pan_db = f"anvio/{PROJECT_NAME}/{PROJECT_NAME}-PAN.db",
        done = touch(f"anvio/{PROJECT_NAME}/pangenome.done")
    benchmark: "logs/anvio/pangenome.benchmark.txt"
    log: "logs/anvio/pangenome.log"
    threads: config.get("pangenome_threads", 32)
    params:
        mcl_inflation = config.get("mcl_inflation", 10),
        min_occurrence = config.get("min_occurrence", 2)
    conda: "envs/anvio.yaml"
    shell:
        "anvi-pan-genome "
        "-g {input.storage} "
        "-n {config[project_name]} "
        "--output-dir anvio/{config[project_name]} "
        "--num-threads {threads} "
        "--use-ncbi-blast "
        "--mcl-inflation {params.mcl_inflation} "
        "--min-occurrence {params.min_occurrence} "
        "2>&1 | tee {log}"

# ======================================================================
# 13. Compute ANI (Average Nucleotide Identity)
# ======================================================================

rule compute_ani:
    """
    Compute ANI between all genomes using PyANI v.0.2.
    """
    input:
        external = "anvio/external-genomes.txt",
        pan_db = f"anvio/{PROJECT_NAME}/{PROJECT_NAME}-PAN.db",
        pangenome_done = f"anvio/{PROJECT_NAME}/pangenome.done"
    output:
        ani_matrix = f"anvio/{PROJECT_NAME}/ANI/ANI_percentage_identity.txt",
        done = touch(f"anvio/{PROJECT_NAME}/ani.done")
    benchmark: "logs/anvio/ani.benchmark.txt"
    log: "logs/anvio/ani.log"
    threads: config.get("ani_threads", 16)
    conda: "envs/anvio.yaml"
    shell:
        "anvi-compute-genome-similarity "
        "-e {input.external} "
        "-p {input.pan_db} "
        "-o anvio/{config[project_name]}/ANI/ "
        "--program pyANI "
        "--num-threads {threads} "
        "2>&1 | tee {log}"

# ======================================================================
# 14. Extract single-copy gene clusters (SCGs)
# ======================================================================

rule extract_scgs:
    """
    Extract single-copy gene clusters: present in ALL genomes,
    exactly ONE gene per genome per cluster.
    Cleans FASTA headers for compatibility with MUSCLE/IQ-TREE.
    """
    input:
        storage = f"anvio/{PROJECT_NAME}_GENOMES.db",
        pan_db = f"anvio/{PROJECT_NAME}/{PROJECT_NAME}-PAN.db",
        pangenome_done = f"anvio/{PROJECT_NAME}/pangenome.done"
    output:
        proteins = f"anvio/{PROJECT_NAME}/concatenated-proteins.fa",
        proteins_clean = f"anvio/{PROJECT_NAME}/concatenated-proteins_clean.fa",
        scg_list = f"anvio/{PROJECT_NAME}/SCGs.txt",
        done = touch(f"anvio/{PROJECT_NAME}/scgs.done")
    benchmark: "logs/anvio/scgs.benchmark.txt"
    log: "logs/anvio/scgs.log"
    conda: "envs/anvio.yaml"
    shell:
        # Get SCGs: present in all genomes, single-copy per genome
        "anvi-get-sequences-for-gene-clusters "
        "-g {input.storage} "
        "-p {input.pan_db} "
        "--min-num-genomes {config[n_genomes]} "
        "--max-num-genomes {config[n_genomes]} "
        "--min-num-genes 1 "
        "--max-num-genes 1 "
        "--concatenate-gene-clusters "
        "-o {output.proteins} "
        "2>&1 | tee {log} && "
        # CRITICAL: Clean FASTA headers for MUSCLE/IQ-TREE compatibility
        # Removes everything after first space in header lines
        "sed -i 's/^\\(>[^[:space:]]*\\)[[:space:]].*/\\1/' {output.proteins} 2>&1 | tee -a {log} && "
        # Create clean copy for safety
        "cp {output.proteins} {output.proteins_clean} && "
        # Save cluster names
        "sqlite3 {input.pan_db} \"SELECT gene_cluster_id FROM gene_clusters "
        "WHERE num_genomes = {config[n_genomes]} AND num_genes = 1;\" "
        "> {output.scg_list} 2>&1 | tee -a {log}"

# ======================================================================
# 15. Phylogenomic tree
# ======================================================================

rule phylogenomic_tree:
    """
    Build phylogenomic tree from concatenated SCG protein sequences.
    Uses MUSCLE v5 for alignment and IQ-TREE for tree building
    with 1000 bootstrap iterations and automatic model selection.
    """
    input:
        proteins = f"anvio/{PROJECT_NAME}/concatenated-proteins_clean.fa",
        scgs_done = f"anvio/{PROJECT_NAME}/scgs.done"
    output:
        tree = f"anvio/{PROJECT_NAME}/phylogenomic-tree.txt",
        done = touch(f"anvio/{PROJECT_NAME}/phylogeny.done")
    benchmark: "logs/anvio/phylogeny.benchmark.txt"
    log: "logs/anvio/phylogeny.log"
    threads: config.get("phylogeny_threads", 32)
    conda: "envs/anvio.yaml"
    shell:
        # Align with MUSCLE v5
        "muscle "
        "-align {input.proteins} "
        "-output anvio/{config[project_name]}/aligned-proteins.fa "
        "-threads {threads} "
        "2>&1 | tee {log} && "
        # Build tree with IQ-TREE (1000 ultrafast bootstrap + SH-aLRT)
        "iqtree "
        "-s anvio/{config[project_name]}/aligned-proteins.fa "
        "-nt {threads} "
        "-bb 1000 "
        "-alrt 1000 "
        "-m MFP "
        "-pre anvio/{config[project_name]}/iqtree "
        "2>&1 | tee -a {log} && "
        "cp anvio/{config[project_name]}/iqtree.treefile {output.tree}"

# ======================================================================
# 16. Functional enrichment analysis
# ======================================================================

rule functional_enrichment:
    """
    Compute functional enrichment in pangenome using GLM
    (Generalized Linear Model) via enrichR.
    """
    input:
        pan_db = f"anvio/{PROJECT_NAME}/{PROJECT_NAME}-PAN.db",
        storage = f"anvio/{PROJECT_NAME}_GENOMES.db",
        pangenome_done = f"anvio/{PROJECT_NAME}/pangenome.done"
    output:
        enrichment = f"anvio/{PROJECT_NAME}/functional_enrichment.txt",
        done = touch(f"anvio/{PROJECT_NAME}/func_enrich.done")
    benchmark: "logs/anvio/func_enrich.benchmark.txt"
    log: "logs/anvio/func_enrich.log"
    threads: config.get("enrichment_threads", 16)
    conda: "envs/anvio.yaml"
    shell:
        "anvi-compute-functional-enrichment-in-pan "
        "-p {input.pan_db} "
        "-g {input.storage} "
        "--annotation-source KOfam "
        "--include-gc-identity-as-function "
        "-o {output.enrichment} "
        "--functional-occurrence-table-output "
        "anvio/{config[project_name]}/functional_occurrence.txt "
        "--num-threads {threads} "
        "2>&1 | tee {log}"

# ======================================================================
# 17. Metabolism estimation
# ======================================================================

rule metabolism_estimation:
    """
    Estimate metabolic capabilities from KEGG annotations.
    Uses module completion threshold to determine pathway presence.
    """
    input:
        external = "anvio/external-genomes.txt",
        kegg_done = expand("anvio/databases/{genome}_kegg.done", genome=GENOME_NAMES)
    output:
        metabolism = f"anvio/{PROJECT_NAME}/metabolism_estimation.txt",
        done = touch(f"anvio/{PROJECT_NAME}/metabolism.done")
    benchmark: "logs/anvio/metabolism.benchmark.txt"
    log: "logs/anvio/metabolism.log"
    threads: config.get("metabolism_threads", 16)
    params:
        threshold = config.get("module_completion_threshold", 0.5)
    conda: "envs/anvio.yaml"
    shell:
        "anvi-estimate-metabolism "
        "-e {input.external} "
        "--module-completion-threshold {params.threshold} "
        "-o {output.metabolism} "
        "--num-threads {threads} "
        "2>&1 | tee {log}"

# ======================================================================
# 18. Summary (export all pangenome data)
# ======================================================================

rule pangenome_summary:
    """
    Generate comprehensive pangenome summary.
    """
    input:
        storage = f"anvio/{PROJECT_NAME}_GENOMES.db",
        pan_db = f"anvio/{PROJECT_NAME}/{PROJECT_NAME}-PAN.db",
        pangenome_done = f"anvio/{PROJECT_NAME}/pangenome.done"
    output:
        done = touch(f"anvio/{PROJECT_NAME}/pangenome_summary.done")
    benchmark: "logs/anvio/summary.benchmark.txt"
    log: "logs/anvio/summary.log"
    conda: "envs/anvio.yaml"
    shell:
        "anvi-summarize "
        "-p {input.pan_db} "
        "-g {input.storage} "
        "-o anvio/{config[project_name]}/summary "
        "-C DEFAULT "
        "2>&1 | tee {log}"
