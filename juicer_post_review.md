# Snakefile_post_review – Regenerate assembly & re-analyze after Juicebox review
# Input: .review.assembly from Juicebox GUI
# Output: Corrected FASTA, new .hic map, QUAST stats, full Hi‑C re-analysis

import os
import glob

configfile: "config_juicer.yaml"

# ======================================================================
# Helpers
# ======================================================================

def find_review_assemblies(wildcards):
    """Find all .review.assembly files for a sample's MAGs."""
    review_files = []
    juicer_dir = f"juicer/{wildcards.sample}"
    for root, dirs, files in os.walk(juicer_dir):
        for f in files:
            if f.endswith('.review.assembly'):
                review_files.append(os.path.join(root, f))
    return review_files

def get_mag_from_review(review_path):
    """Extract MAG name from review assembly path."""
    # juicer/{sample}/{mag_name}/3d-dna/genome.review.assembly -> mag_name
    parts = review_path.split('/')
    if len(parts) >= 3:
        return parts[-3]
    return None

# ======================================================================
# Final targets
# ======================================================================

rule all:
    input:
        expand("post_review/{sample}/all_post_review.done", sample=config["samples"]),

# ======================================================================
# 1. Generate corrected FASTA from .review.assembly
# ======================================================================

rule generate_corrected_fasta:
    """
    Use 3D-DNA's run-asm-pipeline-post-review.sh to generate
    the corrected FASTA file from manual Juicebox edits.
    """
    input:
        review = "juicer/{sample}/{mag_name}/3d-dna/{review_file}",
        original_fasta = "juicer/{sample}/references/{mag_name}/{mag_name}.fa",
        merged_nodups = "juicer/{sample}/{mag_name}/aligned/merged_nodups.txt"
    output:
        corrected_fasta = "post_review/{sample}/{mag_name}/genome_corrected.fasta",
        corrected_assembly = "post_review/{sample}/{mag_name}/genome_corrected.assembly",
        done = touch("post_review/{sample}/{mag_name}/corrected_fasta.done")
    benchmark: "logs/post_review/correct_fasta/{sample}/{mag_name}.benchmark.txt"
    log: "logs/post_review/correct_fasta/{sample}/{mag_name}.log"
    threads: config.get("n_threads", 32)
    params:
        juicer_path = config.get("juicer_path", "auto"),
        out_dir = lambda wc: f"post_review/{wc.sample}/{wc.mag_name}"
    conda: "envs/juicer.yaml"
    shell:
        "if [ \"{params.juicer_path}\" = \"auto\" ]; then "
        "  JUICER_PATH=$(dirname $(dirname $(which juicer.sh))); "
        "else "
        "  JUICER_PATH={params.juicer_path}; "
        "fi && "
        "mkdir -p {params.out_dir} && "
        "cp {input.review} {params.out_dir}/genome.review.assembly && "
        "cd {params.out_dir} && "
        "bash $JUICER_PATH/3d-dna/run-asm-pipeline-post-review.sh "
        "-r genome.review.assembly "
        "$(realpath {input.original_fasta}) "
        "$(realpath {input.merged_nodups}) "
        "2>&1 | tee {log} && "
        # The script outputs: genome.FINAL.fasta, genome.FINAL.assembly
        "if [ -f \"genome.FINAL.fasta\" ]; then "
        "  cp genome.FINAL.fasta {output.corrected_fasta}; "
        "  cp genome.FINAL.assembly {output.corrected_assembly}; "
        "elif [ -f \"genome.fasta\" ]; then "
        "  cp genome.fasta {output.corrected_fasta}; "
        "  cp genome.assembly {output.corrected_assembly}; "
        "fi && "
        "touch {output.done}"

# ======================================================================
# 2. Generate new clean .hic map from corrected FASTA
# ======================================================================

rule generate_corrected_hic:
    """
    Generate a new .hic contact map using the corrected FASTA.
    Follows the same Juicer workflow but with the reviewed assembly.
    """
    input:
        corrected_fasta = "post_review/{sample}/{mag_name}/genome_corrected.fasta",
        merged_nodups = "juicer/{sample}/{mag_name}/aligned/merged_nodups.txt",
        corrected_done = "post_review/{sample}/{mag_name}/corrected_fasta.done"
    output:
        corrected_hic = "post_review/{sample}/{mag_name}/hic/inter_30_corrected.hic",
        done = touch("post_review/{sample}/{mag_name}/corrected_hic.done")
    benchmark: "logs/post_review/correct_hic/{sample}/{mag_name}.benchmark.txt"
    log: "logs/post_review/correct_hic/{sample}/{mag_name}.log"
    threads: config.get("n_threads", 32)
    params:
        juicer_path = config.get("juicer_path", "auto"),
        resolutions = config.get("hic_resolutions", "5000,10000,25000,50000,100000,250000,500000,1000000"),
        mem = config.get("juicer_java_mem", "50g")
    conda: "envs/juicer.yaml"
    shell:
        "if [ \"{params.juicer_path}\" = \"auto\" ]; then "
        "  JUICER_PATH=$(dirname $(dirname $(which juicer.sh))); "
        "else "
        "  JUICER_PATH={params.juicer_path}; "
        "fi && "
        "JTOOLS=\"\" && "
        "if [ -f \"$(which juicer_tools.jar)\" ]; then "
        "  JTOOLS=$(which juicer_tools.jar); "
        "elif [ -f \"$JUICER_PATH/juicer_tools.jar\" ]; then "
        "  JTOOLS=$JUICER_PATH/juicer_tools.jar; "
        "fi && "
        "if [ -z \"$JTOOLS\" ]; then "
        "  echo 'juicer_tools.jar not found' | tee {log}; "
        "  touch {output.done}; exit 0; "
        "fi && "
        "mkdir -p post_review/{wildcards.sample}/{wildcards.mag_name}/hic && "
        # Create chrom.sizes from corrected FASTA
        "samtools faidx {input.corrected_fasta} && "
        "cut -f1,2 {input.corrected_fasta}.fai > post_review/{wildcards.sample}/{wildcards.mag_name}/hic/chrom.sizes && "
        # Generate new .hic
        "java -Xmx{params.mem} -jar $JTOOLS pre "
        "-r {params.resolutions} "
        "{input.merged_nodups} "
        "{output.corrected_hic} "
        "post_review/{wildcards.sample}/{wildcards.mag_name}/hic/chrom.sizes "
        "2>&1 | tee {log} && "
        "touch {output.done}"

# ======================================================================
# 3. QUAST on corrected assembly
# ======================================================================

rule quast_corrected:
    """
    Generate QUAST statistics for the corrected assembly.
    Compare with original if available.
    """
    input:
        corrected_fasta = "post_review/{sample}/{mag_name}/genome_corrected.fasta",
        original_fasta = "juicer/{sample}/references/{mag_name}/{mag_name}.fa",
        corrected_done = "post_review/{sample}/{mag_name}/corrected_fasta.done"
    output:
        html = "post_review/{sample}/{mag_name}/quast/report.html",
        done = touch("post_review/{sample}/{mag_name}/quast.done")
    benchmark: "logs/post_review/quast/{sample}/{mag_name}.benchmark.txt"
    log: "logs/post_review/quast/{sample}/{mag_name}.log"
    threads: config.get("n_threads", 16)
    conda: "envs/quast.yaml"
    params:
        out_dir = lambda wc: f"post_review/{wc.sample}/{wc.mag_name}/quast"
    shell:
        "mkdir -p {params.out_dir} && "
        "quast.py {input.corrected_fasta} {input.original_fasta} "
        "-o {params.out_dir} "
        "--threads {threads} "
        "--fungus "
        "--gene-finding "
        "--rna-finding "
        "--conserved-genes-finding "
        "--circos "
        "--labels corrected,original "
        "--html-report {output.html} "
        "2>&1 | tee {log} && "
        "touch {output.done}"

# ======================================================================
# 4. Assembly statistics (N50, L50, etc.)
# ======================================================================

rule assembly_stats:
    """
    Calculate detailed assembly statistics for the corrected FASTA.
    """
    input:
        corrected_fasta = "post_review/{sample}/{mag_name}/genome_corrected.fasta",
        corrected_done = "post_review/{sample}/{mag_name}/corrected_fasta.done"
    output:
        stats = "post_review/{sample}/{mag_name}/assembly_stats.txt",
        done = touch("post_review/{sample}/{mag_name}/stats.done")
    benchmark: "logs/post_review/stats/{sample}/{mag_name}.benchmark.txt"
    log: "logs/post_review/stats/{sample}/{mag_name}.log"
    conda: "envs/bioinf_v1.yaml"
    shell:
        "python scripts/assembly_stats.py "
        "--fasta {input.corrected_fasta} "
        "--output {output.stats} "
        "--sample {wildcards.mag_name} 2>&1 | tee {log} && "
        "touch {output.done}"

# ======================================================================
# 5. Re-run full Hi‑C analysis on corrected map
# ======================================================================

if config.get("run_hic_analysis", True):
    
    rule post_review_compartments:
        input:
            hic = "post_review/{sample}/{mag_name}/hic/inter_30_corrected.hic",
            hic_done = "post_review/{sample}/{mag_name}/corrected_hic.done"
        output:
            eigen = "post_review/{sample}/{mag_name}/analysis/compartments/eigenvectors.txt",
            done = touch("post_review/{sample}/{mag_name}/analysis/compartments.done")
        benchmark: "logs/post_review/compartments/{sample}/{mag_name}.benchmark.txt"
        log: "logs/post_review/compartments/{sample}/{mag_name}.log"
        params:
            juicer_path = config.get("juicer_path", "auto"),
            mem = config.get("analysis_mem", config.get("juicer_java_mem", "30g"))
        conda: "envs/juicer.yaml"
        shell:
            "if [ \"{params.juicer_path}\" = \"auto\" ]; then "
            "  JUICER_PATH=$(dirname $(dirname $(which juicer.sh))); "
            "else "
            "  JUICER_PATH={params.juicer_path}; "
            "fi && "
            "JTOOLS=\"\" && "
            "if [ -f \"$(which juicer_tools.jar)\" ]; then "
            "  JTOOLS=$(which juicer_tools.jar); "
            "elif [ -f \"$JUICER_PATH/juicer_tools.jar\" ]; then "
            "  JTOOLS=$JUICER_PATH/juicer_tools.jar; "
            "fi && "
            "if [ -z \"$JTOOLS\" ]; then "
            "  touch {output.done}; exit 0; "
            "fi && "
            "mkdir -p post_review/{wildcards.sample}/{wildcards.mag_name}/analysis/compartments && "
            "java -Xmx{params.mem} -jar $JTOOLS eigenvector "
            "-p KR "
            "{input.hic} "
            "post_review/{wildcards.sample}/{wildcards.mag_name}/analysis/compartments/eigenvectors "
            "2>&1 | tee {log} && "
            "if [ -f \"post_review/{wildcards.sample}/{wildcards.mag_name}/analysis/compartments/eigenvectors_KR.txt\" ]; then "
            "  mv post_review/{wildcards.sample}/{wildcards.mag_name}/analysis/compartments/eigenvectors_KR.txt {output.eigen}; "
            "fi && "
            "touch {output.done}"

    rule post_review_distance_plot:
        input:
            hic = "post_review/{sample}/{mag_name}/hic/inter_30_corrected.hic",
            hic_done = "post_review/{sample}/{mag_name}/corrected_hic.done"
        output:
            plot = "post_review/{sample}/{mag_name}/analysis/qc/distance_plot_corrected.png",
            done = touch("post_review/{sample}/{mag_name}/analysis/qc.done")
        benchmark: "logs/post_review/distance_plot/{sample}/{mag_name}.benchmark.txt"
        log: "logs/post_review/distance_plot/{sample}/{mag_name}.log"
        params:
            juicer_path = config.get("juicer_path", "auto"),
            mem = config.get("analysis_mem", config.get("juicer_java_mem", "30g"))
        conda: "envs/juicer.yaml"
        shell:
            "if [ \"{params.juicer_path}\" = \"auto\" ]; then "
            "  JUICER_PATH=$(dirname $(dirname $(which juicer.sh))); "
            "else "
            "  JUICER_PATH={params.juicer_path}; "
            "fi && "
            "JTOOLS=\"\" && "
            "if [ -f \"$(which juicer_tools.jar)\" ]; then "
            "  JTOOLS=$(which juicer_tools.jar); "
            "elif [ -f \"$JUICER_PATH/juicer_tools.jar\" ]; then "
            "  JTOOLS=$JUICER_PATH/juicer_tools.jar; "
            "fi && "
            "if [ -z \"$JTOOLS\" ]; then "
            "  touch {output.done}; exit 0; "
            "fi && "
            "mkdir -p post_review/{wildcards.sample}/{wildcards.mag_name}/analysis/qc && "
            "python scripts/distance_plot.py "
            "--hic {input.hic} "
            "--output {output.plot} "
            "--juicer-tools $JTOOLS "
            "--sample {wildcards.mag_name}_corrected 2>&1 | tee {log} && "
            "touch {output.done}"

    # Additional analysis steps (domains, loops, APA) follow same pattern
    # but point to post_review/{sample}/{mag_name}/analysis/ instead of hic_analysis/

    rule post_review_analysis_complete:
        input:
            compartments = "post_review/{sample}/{mag_name}/analysis/compartments.done",
            qc = "post_review/{sample}/{mag_name}/analysis/qc.done"
        output:
            done = touch("post_review/{sample}/{mag_name}/analysis/analysis_complete.done")
        shell: "touch {output.done}"

# ======================================================================
# 6. Comparison report: original vs corrected
# ======================================================================

rule comparison_report:
    """
    Generate a comparison report between original and corrected assemblies.
    """
    input:
        original_stats = "juicer/{sample}/{mag_name}/aligned/statistics.txt",
        original_quast = "results/{sample}/quast/report.html" if config.get("run_quast", True) else [],
        corrected_stats = "post_review/{sample}/{mag_name}/assembly_stats.txt",
        corrected_quast = "post_review/{sample}/{mag_name}/quast/report.html",
        corrected_done = "post_review/{sample}/{mag_name}/corrected_fasta.done",
        stats_done = "post_review/{sample}/{mag_name}/stats.done",
        quast_done = "post_review/{sample}/{mag_name}/quast.done"
    output:
        html = "post_review/{sample}/{mag_name}/comparison_report.html",
        done = touch("post_review/{sample}/{mag_name}/comparison.done")
    benchmark: "logs/post_review/comparison/{sample}/{mag_name}.benchmark.txt"
    log: "logs/post_review/comparison/{sample}/{mag_name}.log"
    conda: "envs/bioinf_v1.yaml"
    shell:
        "python scripts/generate_comparison_report.py "
        "--original-stats {input.original_stats} "
        "--corrected-stats {input.corrected_stats} "
        "--corrected-quast {input.corrected_quast} "
        "--output {output.html} "
        "--sample {wildcards.mag_name} 2>&1 | tee {log} && "
        "touch {output.done}"

# ======================================================================
# 7. Per-MAG completion
# ======================================================================

rule mag_post_review_complete:
    input:
        corrected_fasta = "post_review/{sample}/{mag_name}/corrected_fasta.done",
        corrected_hic = "post_review/{sample}/{mag_name}/corrected_hic.done",
        quast = "post_review/{sample}/{mag_name}/quast.done",
        stats = "post_review/{sample}/{mag_name}/stats.done",
        analysis = "post_review/{sample}/{mag_name}/analysis/analysis_complete.done" if config.get("run_hic_analysis", True) else [],
        comparison = "post_review/{sample}/{mag_name}/comparison.done"
    output:
        done = touch("post_review/{sample}/{mag_name}/post_review_complete.done")
    shell: "touch {output.done}"

# ======================================================================
# 8. Sample-level completion
# ======================================================================

def get_mag_names_for_sample(sample):
    """Get MAG names from Juicer output directory."""
    mag_list = f"juicer/{sample}/mag_list.txt"
    if os.path.exists(mag_list):
        with open(mag_list) as f:
            return [line.strip() for line in f if line.strip()]
    return []

rule sample_post_review_complete:
    input:
        lambda wc: expand(
            "post_review/{sample}/{mag_name}/post_review_complete.done",
            sample=wc.sample,
            mag_name=get_mag_names_for_sample(wc.sample)
        )
    output:
        done = touch("post_review/{sample}/all_post_review.done")
    shell: "touch {output.done}"
