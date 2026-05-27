# Snakefile_hic_analysis – Downstream Hi‑C Analysis for Fungal MAGs
# Takes .hic files from Juicer pipeline, performs normalization,
# feature annotation, compartment analysis, and QC.

import os

configfile: "config_hic_analysis.yaml"

# ======================================================================
# Helper: get all .hic files for a sample
# ======================================================================

def get_hic_files(wildcards):
    """Find all .hic files for a given sample."""
    hic_dir = f"juicer/{wildcards.sample}"
    hic_files = []
    for root, dirs, files in os.walk(hic_dir):
        for f in files:
            if f.endswith('.hic'):
                hic_files.append(os.path.join(root, f))
    return hic_files

def get_mag_name(hic_path):
    """Extract MAG name from .hic file path."""
    # juicer/{sample}/{mag_name}/aligned/inter_30.hic -> mag_name
    parts = hic_path.split('/')
    if len(parts) >= 3:
        return parts[-3]
    return os.path.basename(os.path.dirname(os.path.dirname(hic_path)))

# ======================================================================
# Final targets
# ======================================================================

rule all:
    input:
        expand("hic_analysis/{sample}/all_analysis.done", sample=config["samples"]),

# ======================================================================
# 1. Prepare analysis directory per MAG
# ======================================================================

rule prepare_analysis_dir:
    """
    Create analysis directory structure for each MAG's .hic file.
    """
    input:
        hic = "juicer/{sample}/{mag_name}/aligned/inter_30.hic"
    output:
        done = touch("hic_analysis/{sample}/{mag_name}/prepared.done")
    benchmark: "logs/hic_analysis/prepare/{sample}/{mag_name}.benchmark.txt"
    log: "logs/hic_analysis/prepare/{sample}/{mag_name}.log"
    run:
        analysis_dir = f"hic_analysis/{wildcards.sample}/{wildcards.mag_name}"
        for subdir in ['compartments', 'domains', 'loops', 'qc', 'matrices', 'apa']:
            os.makedirs(os.path.join(analysis_dir, subdir), exist_ok=True)

# ======================================================================
# 2. Apply multiple normalizations to .hic
# ======================================================================

rule apply_normalizations:
    """
    Generate normalized .hic files using multiple methods:
    - KR (Knight-Ruiz) – iterative matrix balancing
    - VC (Vanilla Coverage) – coverage normalization
    - VC_SQRT – square root of VC
    - SCALE – multiplicative normalization

    Uses juicer_tools to create normalized versions.
    """
    input:
        hic = "juicer/{sample}/{mag_name}/aligned/inter_30.hic",
        prepared = "hic_analysis/{sample}/{mag_name}/prepared.done"
    output:
        kr = "hic_analysis/{sample}/{mag_name}/normalized/kr.hic",
        vc = "hic_analysis/{sample}/{mag_name}/normalized/vc.hic",
        vc_sqrt = "hic_analysis/{sample}/{mag_name}/normalized/vc_sqrt.hic",
        done = touch("hic_analysis/{sample}/{mag_name}/normalizations.done")
    benchmark: "logs/hic_analysis/normalize/{sample}/{mag_name}.benchmark.txt"
    log: "logs/hic_analysis/normalize/{sample}/{mag_name}.log"
    threads: config.get("threads", 16)
    params:
        juicer_path = config.get("juicer_path", "auto"),
        resolutions = config.get("hic_resolutions", "5000,10000,25000,50000,100000,250000,500000,1000000"),
        mem = config.get("juicer_mem", "30g")
    conda: "envs/juicer.yaml"
    shell:
        # Find juicer_tools.jar
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
        "  touch {output.done}; "
        "  exit 0; "
        "fi && "
        # Create output directory
        "mkdir -p hic_analysis/{wildcards.sample}/{wildcards.mag_name}/normalized && "
        # KR normalization
        "echo 'Running KR normalization...' | tee {log} && "
        "java -Xmx{params.mem} -jar $JTOOLS addNorm -r {params.resolutions} "
        "-N KR "
        "{input.hic} "
        "{output.kr} 2>&1 | tee -a {log} && "
        # VC normalization
        "echo 'Running VC normalization...' | tee -a {log} && "
        "java -Xmx{params.mem} -jar $JTOOLS addNorm -r {params.resolutions} "
        "-N VC "
        "{input.hic} "
        "{output.vc} 2>&1 | tee -a {log} && "
        # VC_SQRT normalization
        "echo 'Running VC_SQRT normalization...' | tee -a {log} && "
        "java -Xmx{params.mem} -jar $JTOOLS addNorm -r {params.resolutions} "
        "-N VC_SQRT "
        "{input.hic} "
        "{output.vc_sqrt} 2>&1 | tee -a {log} && "
        "touch {output.done}"

# ======================================================================
# 3. Compartment analysis (A/B compartments)
# ======================================================================

rule eigenvector_compartments:
    """
    Calculate A/B compartments at multiple resolutions using eigenvector.
    A compartments = open/active chromatin
    B compartments = closed/inactive chromatin
    """
    input:
        hic = "hic_analysis/{sample}/{mag_name}/normalized/kr.hic",
        normalizations = "hic_analysis/{sample}/{mag_name}/normalizations.done"
    output:
        eigen = "hic_analysis/{sample}/{mag_name}/compartments/eigenvectors.txt",
        done = touch("hic_analysis/{sample}/{mag_name}/compartments.done")
    benchmark: "logs/hic_analysis/compartments/{sample}/{mag_name}.benchmark.txt"
    log: "logs/hic_analysis/compartments/{sample}/{mag_name}.log"
    params:
        juicer_path = config.get("juicer_path", "auto"),
        resolutions = "10000,25000,50000,100000",
        mem = config.get("juicer_mem", "30g")
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
        "mkdir -p hic_analysis/{wildcards.sample}/{wildcards.mag_name}/compartments && "
        "echo 'Calculating eigenvectors...' | tee {log} && "
        "java -Xmx{params.mem} -jar $JTOOLS eigenvector "
        "-p KR "
        "{input.hic} "
        "hic_analysis/{wildcards.sample}/{wildcards.mag_name}/compartments/eigenvectors "
        "2>&1 | tee -a {log} && "
        # Rename output to match expected filename
        "if [ -f \"hic_analysis/{wildcards.sample}/{wildcards.mag_name}/compartments/eigenvectors_KR.txt\" ]; then "
        "  mv hic_analysis/{wildcards.sample}/{wildcards.mag_name}/compartments/eigenvectors_KR.txt {output.eigen}; "
        "fi && "
        "touch {output.done}"

# ======================================================================
# 4. Arrowhead – contact domain detection
# ======================================================================

rule arrowhead_domains:
    """
    Detect contact domains (TADs in mammals, domain-like structures in fungi)
    using Juicer's Arrowhead algorithm.
    """
    input:
        hic = "hic_analysis/{sample}/{mag_name}/normalized/kr.hic",
        normalizations = "hic_analysis/{sample}/{mag_name}/normalizations.done"
    output:
        domains = "hic_analysis/{sample}/{mag_name}/domains/arrowhead_domains.bedpe",
        done = touch("hic_analysis/{sample}/{mag_name}/domains.done")
    benchmark: "logs/hic_analysis/domains/{sample}/{mag_name}.benchmark.txt"
    log: "logs/hic_analysis/domains/{sample}/{mag_name}.log"
    params:
        juicer_path = config.get("juicer_path", "auto"),
        resolutions = "5000,10000,25000",
        mem = config.get("juicer_mem", "30g")
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
        "mkdir -p hic_analysis/{wildcards.sample}/{wildcards.mag_name}/domains && "
        "echo 'Running Arrowhead domain detection...' | tee {log} && "
        "java -Xmx{params.mem} -jar $JTOOLS arrowhead "
        "-r {params.resolutions} "
        "-k KR "
        "{input.hic} "
        "hic_analysis/{wildcards.sample}/{wildcards.mag_name}/domains/arrowhead "
        "2>&1 | tee -a {log} && "
        # Merge resolutions into one BEDPE
        "python scripts/merge_arrowhead_domains.py "
        "hic_analysis/{wildcards.sample}/{wildcards.mag_name}/domains/arrowhead "
        "{output.domains} 2>&1 | tee -a {log} && "
        "touch {output.done}"

# ======================================================================
# 5. HiCCUPS – loop detection
# ======================================================================

rule hiccups_loops:
    """
    Detect chromatin loops using Juicer's HiCCUPS algorithm.
    Loops indicate long-range interactions.
    """
    input:
        hic = "hic_analysis/{sample}/{mag_name}/normalized/kr.hic",
        normalizations = "hic_analysis/{sample}/{mag_name}/normalizations.done"
    output:
        loops = "hic_analysis/{sample}/{mag_name}/loops/hiccups_loops.bedpe",
        done = touch("hic_analysis/{sample}/{mag_name}/loops.done")
    benchmark: "logs/hic_analysis/loops/{sample}/{mag_name}.benchmark.txt"
    log: "logs/hic_analysis/loops/{sample}/{mag_name}.log"
    params:
        juicer_path = config.get("juicer_path", "auto"),
        resolutions = "5000,10000,25000",
        mem = config.get("juicer_mem", "30g")
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
        "mkdir -p hic_analysis/{wildcards.sample}/{wildcards.mag_name}/loops && "
        "echo 'Running HiCCUPS loop detection...' | tee {log} && "
        "java -Xmx{params.mem} -jar $JTOOLS hiccups "
        "-r {params.resolutions} "
        "-k KR "
        "--cpu --threads {threads} "
        "{input.hic} "
        "hic_analysis/{wildcards.sample}/{wildcards.mag_name}/loops/hiccups "
        "2>&1 | tee -a {log} && "
        "touch {output.done}"

# ======================================================================
# 6. APA – Aggregate Peak Analysis
# ======================================================================

rule apa_analysis:
    """
    Perform Aggregate Peak Analysis on detected loops.
    APA validates loop calls by aggregating signal across all loops.
    """
    input:
        hic = "hic_analysis/{sample}/{mag_name}/normalized/kr.hic",
        loops = "hic_analysis/{sample}/{mag_name}/loops/hiccups_loops.bedpe",
        loops_done = "hic_analysis/{sample}/{mag_name}/loops.done",
        normalizations = "hic_analysis/{sample}/{mag_name}/normalizations.done"
    output:
        apa = "hic_analysis/{sample}/{mag_name}/apa/apa_results.txt",
        plots = "hic_analysis/{sample}/{mag_name}/apa/apa_plots.pdf",
        done = touch("hic_analysis/{sample}/{mag_name}/apa.done")
    benchmark: "logs/hic_analysis/apa/{sample}/{mag_name}.benchmark.txt"
    log: "logs/hic_analysis/apa/{sample}/{mag_name}.log"
    params:
        juicer_path = config.get("juicer_path", "auto"),
        resolutions = "5000,10000",
        mem = config.get("juicer_mem", "30g")
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
        "if [ -z \"$JTOOLS\" ] || [ ! -s {input.loops} ]; then "
        "  touch {output.done}; exit 0; "
        "fi && "
        "mkdir -p hic_analysis/{wildcards.sample}/{wildcards.mag_name}/apa && "
        "echo 'Running APA analysis...' | tee {log} && "
        "java -Xmx{params.mem} -jar $JTOOLS apa "
        "-r {params.resolutions} "
        "-k KR "
        "{input.hic} "
        "{input.loops} "
        "hic_analysis/{wildcards.sample}/{wildcards.mag_name}/apa/apa "
        "2>&1 | tee -a {log} && "
        "touch {output.done}"

# ======================================================================
# 7. Distance vs. Contact Frequency Plot (QC)
# ======================================================================

rule distance_plot:
    """
    Generate distance-dependent contact frequency plot.
    This is a fundamental QC metric for Hi‑C data quality.
    """
    input:
        hic = "juicer/{sample}/{mag_name}/aligned/inter_30.hic",
        prepared = "hic_analysis/{sample}/{mag_name}/prepared.done"
    output:
        plot = "hic_analysis/{sample}/{mag_name}/qc/distance_plot.png",
        done = touch("hic_analysis/{sample}/{mag_name}/distance_plot.done")
    benchmark: "logs/hic_analysis/distance_plot/{sample}/{mag_name}.benchmark.txt"
    log: "logs/hic_analysis/distance_plot/{sample}/{mag_name}.log"
    conda: "envs/juicer.yaml"
    params:
        juicer_path = config.get("juicer_path", "auto"),
        mem = config.get("juicer_mem", "30g")
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
        "mkdir -p hic_analysis/{wildcards.sample}/{wildcards.mag_name}/qc && "
        "echo 'Generating distance plot...' | tee {log} && "
        "python scripts/distance_plot.py "
        "--hic {input.hic} "
        "--output {output.plot} "
        "--juicer-tools $JTOOLS "
        "--sample {wildcards.mag_name} 2>&1 | tee -a {log} && "
        "touch {output.done}"

# ======================================================================
# 8. Extract interaction matrices
# ======================================================================

rule extract_matrices:
    """
    Extract raw and normalized interaction matrices at key resolutions.
    Useful for custom analysis and visualization.
    """
    input:
        hic = "hic_analysis/{sample}/{mag_name}/normalized/kr.hic",
        normalizations = "hic_analysis/{sample}/{mag_name}/normalizations.done"
    output:
        done = touch("hic_analysis/{sample}/{mag_name}/matrices.done")
    benchmark: "logs/hic_analysis/matrices/{sample}/{mag_name}.benchmark.txt"
    log: "logs/hic_analysis/matrices/{sample}/{mag_name}.log"
    params:
        juicer_path = config.get("juicer_path", "auto"),
        resolutions = "10000,25000,50000",
        mem = config.get("juicer_mem", "30g")
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
        "MATRIX_DIR=hic_analysis/{wildcards.sample}/{wildcards.mag_name}/matrices && "
        "mkdir -p $MATRIX_DIR && "
        "echo 'Extracting interaction matrices...' | tee {log} && "
        "for RES in {params.resolutions//,/ }; do "
        "  java -Xmx{params.mem} -jar $JTOOLS dump observed KR "
        "  {input.hic} chr1 chr1 BP $RES > $MATRIX_DIR/chr1_chr1_$RES.matrix 2>&1 | tee -a {log}; "
        "done && "
        "touch {output.done}"

# ======================================================================
# 9. Completion marker per MAG
# ======================================================================

rule mag_analysis_complete:
    input:
        compartments = "hic_analysis/{sample}/{mag_name}/compartments.done",
        domains = "hic_analysis/{sample}/{mag_name}/domains.done",
        loops = "hic_analysis/{sample}/{mag_name}/loops.done",
        apa = "hic_analysis/{sample}/{mag_name}/apa.done",
        distance = "hic_analysis/{sample}/{mag_name}/distance_plot.done",
        matrices = "hic_analysis/{sample}/{mag_name}/matrices.done"
    output:
        done = touch("hic_analysis/{sample}/{mag_name}/analysis_complete.done")
    shell: "touch {output.done}"

# ======================================================================
# 10. Sample-level completion
# ======================================================================

rule sample_analysis_complete:
    input:
        lambda wc: expand(
            "hic_analysis/{sample}/{mag_name}/analysis_complete.done",
            sample=wc.sample,
            mag_name=get_mag_names_for_sample(wc.sample)
        )
    output:
        done = touch("hic_analysis/{sample}/all_analysis.done")
    shell: "touch {output.done}"

def get_mag_names_for_sample(sample):
    """Get MAG names from the Juicer output directory."""
    mag_list = f"juicer/{sample}/mag_list.txt"
    if os.path.exists(mag_list):
        with open(mag_list) as f:
            return [line.strip() for line in f if line.strip()]
    return []
