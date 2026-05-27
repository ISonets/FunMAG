# Snakefile – Enhanced Juicer Hi‑C Pipeline
# v4 – Cleaned: no blacklist, optional phasing, -S early mode

import os

configfile: "config_juicer.yaml"

# ======================================================================
# Final targets
# ======================================================================

rule all:
    input:
        expand("juicer/{sample}/all_hic_maps.done", sample=config["samples"]),
        expand("juicer/{sample}/library_complexity_summary.tsv", sample=config["samples"]),

# ======================================================================
# 1. Prepare Hi‑C reads
# ======================================================================

rule prepare_hic_reads:
    input:
        r1 = config["data_dir"] + "/{sample}_HiC_R1.fastq.gz",
        r2 = config["data_dir"] + "/{sample}_HiC_R2.fastq.gz"
    output:
        r1 = "juicer/{sample}/fastq/Hi_C_R1.fastq.gz",
        r2 = "juicer/{sample}/fastq/Hi_C_R2.fastq.gz"
    benchmark: "logs/juicer/prepare_hic/{sample}.benchmark.txt"
    log: "logs/juicer/prepare_hic/{sample}.log"
    shell:
        "mkdir -p juicer/{wildcards.sample}/fastq && "
        "if [ ! -f {output.r1} ]; then ln -sf $(realpath {input.r1}) {output.r1}; fi && "
        "if [ ! -f {output.r2} ]; then ln -sf $(realpath {input.r2}) {output.r2}; fi 2>&1 | tee {log}"

# ======================================================================
# 2. Collect fungal MAGs
# ======================================================================

rule collect_mags_for_juicer:
    input: "results/{sample}/fungal_bins_final.lst"
    output:
        done = touch("juicer/{sample}/mags_collected.done"),
        mag_list = "juicer/{sample}/mag_list.txt"
    benchmark: "logs/juicer/collect_mags/{sample}.benchmark.txt"
    log: "logs/juicer/collect_mags/{sample}.log"
    run:
        mag_dir = f"juicer/{wildcards.sample}/fungal_mags"
        os.makedirs(mag_dir, exist_ok=True)

        mags = []
        if os.path.exists(input[0]):
            with open(input[0]) as f:
                for line in f:
                    src = line.strip()
                    if src and os.path.exists(src):
                        mag_name = os.path.basename(src)
                        mag_name = mag_name.replace('.fa', '').replace('.fna', '').replace('.fasta', '')
                        dest = os.path.join(mag_dir, f"{mag_name}.fa")
                        if not os.path.exists(dest):
                            try:
                                os.link(src, dest)
                            except OSError:
                                import shutil
                                shutil.copy2(src, dest)
                        mags.append(mag_name)

        with open(output.mag_list, 'w') as f:
            for m in mags:
                f.write(m + '\n')

# ======================================================================
# 3. Prepare reference for each MAG
# ======================================================================

rule prepare_reference:
    input: "juicer/{sample}/fungal_mags/{mag_name}.fa"
    output:
        index_done = touch("juicer/{sample}/references/{mag_name}/index.done"),
        sites = "juicer/{sample}/references/{mag_name}/{mag_name}_{enzyme}.txt"
    params:
        enzyme = config["restriction_enzyme"],
        juicer_path = config.get("juicer_path", "auto")
    benchmark: "logs/juicer/prepare_ref/{sample}/{mag_name}.benchmark.txt"
    log: "logs/juicer/prepare_ref/{sample}/{mag_name}.log"
    conda: "envs/juicer.yaml"
    shell:
        "REF_DIR=juicer/{wildcards.sample}/references/{wildcards.mag_name} && "
        "mkdir -p $REF_DIR && "
        "cp {input} $REF_DIR/{wildcards.mag_name}.fa && "
        "bwa index $REF_DIR/{wildcards.mag_name}.fa 2>&1 | tee {log} && "
        "if [ \"{params.juicer_path}\" = \"auto\" ]; then "
        "  JUICER_PATH=$(dirname $(dirname $(which juicer.sh))); "
        "else "
        "  JUICER_PATH={params.juicer_path}; "
        "fi && "
        "python $JUICER_PATH/misc/generate_site_positions.py "
        "{params.enzyme} {wildcards.mag_name} "
        "$REF_DIR/{wildcards.mag_name}.fa 2>&1 | tee -a {log} && "
        "if [ -f \"{wildcards.mag_name}_{params.enzyme}.txt\" ]; then "
        "  mv {wildcards.mag_name}_{params.enzyme}.txt {output.sites}; "
        "elif [ -f \"$REF_DIR/{wildcards.mag_name}_{params.enzyme}.txt\" ]; then "
        "  cp $REF_DIR/{wildcards.mag_name}_{params.enzyme}.txt {output.sites}; "
        "fi && "
        "touch {output.index_done}"

# ======================================================================
# 4. Run Juicer (-S early mode)
# ======================================================================

rule run_juicer:
    """
    Run Juicer in -S early mode (stop after merged_nodups.txt).
    This preserves the deduplicated contact file without auto-generating .hic.
    We manually generate .hic at custom resolutions afterward.
    """
    input:
        mag = "juicer/{sample}/references/{mag_name}/{mag_name}.fa",
        index_done = "juicer/{sample}/references/{mag_name}/index.done",
        sites = "juicer/{sample}/references/{mag_name}/{mag_name}_{enzyme}.txt",
        r1 = "juicer/{sample}/fastq/Hi_C_R1.fastq.gz",
        r2 = "juicer/{sample}/fastq/Hi_C_R2.fastq.gz"
    output:
        merged_nodups = "juicer/{sample}/{mag_name}/aligned/merged_nodups.txt",
        stats = "juicer/{sample}/{mag_name}/aligned/statistics.txt",
        done = touch("juicer/{sample}/{mag_name}/juicer.done")
    benchmark: "logs/juicer/run/{sample}/{mag_name}.benchmark.txt"
    log: "logs/juicer/run/{sample}/{mag_name}.log"
    threads: config["n_threads"]
    params:
        enzyme = config["restriction_enzyme"],
        juicer_path = config.get("juicer_path", "auto")
    conda: "envs/juicer.yaml"
    shell:
        "if [ \"{params.juicer_path}\" = \"auto\" ]; then "
        "  export JUICER_PATH=$(dirname $(dirname $(which juicer.sh))); "
        "else "
        "  export JUICER_PATH={params.juicer_path}; "
        "fi && "
        "OUT_DIR=juicer/{wildcards.sample}/{wildcards.mag_name} && "
        "mkdir -p $OUT_DIR && "
        "cd $OUT_DIR && "
        "REF_ABS=$(realpath ../references/{wildcards.mag_name}/{wildcards.mag_name}.fa) && "
        "SITES_ABS=$(realpath ../references/{wildcards.mag_name}/{wildcards.mag_name}_{params.enzyme}.txt) && "
        # Run Juicer in -S early mode
        "bash $JUICER_PATH/CPU/juicer.sh "
        "-S early "
        "-z $REF_ABS "
        "-g {wildcards.mag_name} "
        "-s {params.enzyme} "
        "-y $SITES_ABS "
        "-D $JUICER_PATH/CPU "
        "-p $(realpath .) "
        "-t {threads} "
        "2>&1 | tee {log} && "
        # Statistics
        "if [ -f \"$JUICER_PATH/CPU/scripts/common/statistics.pl\" ]; then "
        "  perl $JUICER_PATH/CPU/scripts/common/statistics.pl "
        "  -s $SITES_ABS "
        "  aligned/merged_nodups.txt > aligned/statistics.txt 2>&1 | tee -a {log}; "
        "fi && "
        "touch {output.done}"

# ======================================================================
# 5. Generate .hic file at custom resolutions
# ======================================================================

rule generate_hic:
    """
    Generate .hic contact map at custom resolutions appropriate for fungal genomes.
    Runs AFTER Juicer -S early has produced merged_nodups.txt.
    """
    input:
        merged_nodups = "juicer/{sample}/{mag_name}/aligned/merged_nodups.txt",
        juicer_done = "juicer/{sample}/{mag_name}/juicer.done",
        mag = "juicer/{sample}/references/{mag_name}/{mag_name}.fa"
    output:
        hic = "juicer/{sample}/{mag_name}/aligned/inter_30.hic"
    benchmark: "logs/juicer/generate_hic/{sample}/{mag_name}.benchmark.txt"
    log: "logs/juicer/generate_hic/{sample}/{mag_name}.log"
    params:
        resolutions = config.get("hic_resolutions", "5000,10000,25000,50000,100000,250000,500000,1000000"),
        juicer_mem = config.get("juicer_java_mem", "50g"),
        juicer_path = config.get("juicer_path", "auto")
    conda: "envs/juicer.yaml"
    shell:
        "if [ \"{params.juicer_path}\" = \"auto\" ]; then "
        "  export JUICER_PATH=$(dirname $(dirname $(which juicer.sh))); "
        "else "
        "  export JUICER_PATH={params.juicer_path}; "
        "fi && "
        # Find juicer_tools.jar
        "if [ -f \"$(which juicer_tools.jar)\" ]; then "
        "  JTOOLS=$(which juicer_tools.jar); "
        "elif [ -f \"$JUICER_PATH/juicer_tools.jar\" ]; then "
        "  JTOOLS=$JUICER_PATH/juicer_tools.jar; "
        "else "
        "  JTOOLS=''; "
        "fi && "
        "if [ -z \"$JTOOLS\" ]; then "
        "  echo 'juicer_tools.jar not found, skipping .hic generation' | tee {log}; "
        "  touch {output.hic}; "
        "else "
        "  REF_ABS=$(realpath {input.mag}) && "
        "  java -Xmx{params.juicer_mem} -jar $JTOOLS pre "
        "  -r {params.resolutions} "
        "  {input.merged_nodups} "
        "  {output.hic} "
        "  $REF_ABS.fai 2>&1 | tee {log}; "
        "fi"

# ======================================================================
# 6. Library complexity QC
# ======================================================================

rule library_complexity:
    input:
        merged_nodups = "juicer/{sample}/{mag_name}/aligned/merged_nodups.txt",
        juicer_done = "juicer/{sample}/{mag_name}/juicer.done"
    output:
        metrics = "juicer/{sample}/{mag_name}/qc/library_complexity.txt"
    benchmark: "logs/juicer/complexity/{sample}/{mag_name}.benchmark.txt"
    log: "logs/juicer/complexity/{sample}/{mag_name}.log"
    threads: 4
    conda: "envs/juicer.yaml"
    params:
        out_dir = lambda wc: f"juicer/{wc.sample}/{wc.mag_name}/qc"
    shell:
        "mkdir -p {params.out_dir} && "
        "python scripts/calculate_library_complexity.py "
        "--merged-nodups {input.merged_nodups} "
        "--output {output.metrics} "
        "--sample {wildcards.mag_name} 2>&1 | tee {log}"

def get_mags_from_sample(sample):
    """Helper to get MAG names for a sample."""
    mag_list = f"juicer/{sample}/mag_list.txt"
    if os.path.exists(mag_list):
        with open(mag_list) as f:
            return [line.strip() for line in f if line.strip()]
    return []

rule aggregate_complexity:
    input:
        lambda wc: expand(
            "juicer/{sample}/{mag_name}/qc/library_complexity.txt",
            sample=wc.sample,
            mag_name=get_mags_from_sample(wc.sample)
        )
    output:
        summary = "juicer/{sample}/library_complexity_summary.tsv"
    benchmark: "logs/juicer/complexity_summary/{sample}.benchmark.txt"
    log: "logs/juicer/complexity_summary/{sample}.log"
    run:
        header_written = False
        with open(output.summary, 'w') as out:
            for f in input:
                if os.path.exists(f):
                    with open(f) as inf:
                        lines = inf.readlines()
                        if lines:
                            if not header_written:
                                out.write(lines[0])
                                header_written = True
                            out.writelines(lines[1:])

# ======================================================================
# 7. Optional: Diploid phasing (default: off)
# ======================================================================

if config.get("run_phasing", False):
    rule extract_snps:
        input:
            mag = "juicer/{sample}/references/{mag_name}/{mag_name}.fa",
            index_done = "juicer/{sample}/references/{mag_name}/index.done"
        output:
            vcf = "juicer/{sample}/{mag_name}/phasing/snps.vcf"
        benchmark: "logs/juicer/snp_calling/{sample}/{mag_name}.benchmark.txt"
        log: "logs/juicer/snp_calling/{sample}/{mag_name}.log"
        threads: 8
        conda: "envs/juicer.yaml"
        params:
            out_dir = lambda wc: f"juicer/{wc.sample}/{wc.mag_name}/phasing"
        shell:
            "mkdir -p {params.out_dir} && "
            "freebayes -f {input.mag} -p 1 --pooled-continuous "
            "{input.mag} > {output.vcf} 2>&1 | tee {log} || "
            "touch {output.vcf}"

    rule phase_diploid:
        input:
            vcf = "juicer/{sample}/{mag_name}/phasing/snps.vcf",
            merged_nodups = "juicer/{sample}/{mag_name}/aligned/merged_nodups.txt",
            mag = "juicer/{sample}/references/{mag_name}/{mag_name}.fa",
            juicer_done = "juicer/{sample}/{mag_name}/juicer.done"
        output:
            haplotype1 = "juicer/{sample}/{mag_name}/phasing/haplotype_1.fa",
            haplotype2 = "juicer/{sample}/{mag_name}/phasing/haplotype_2.fa",
            phased_vcf = "juicer/{sample}/{mag_name}/phasing/phased.vcf"
        benchmark: "logs/juicer/phasing/{sample}/{mag_name}.benchmark.txt"
        log: "logs/juicer/phasing/{sample}/{mag_name}.log"
        threads: config["n_threads"]
        conda: "envs/juicer.yaml"
        params:
            out_dir = lambda wc: f"juicer/{wc.sample}/{wc.mag_name}/phasing"
        shell:
            "mkdir -p {params.out_dir} && "
            "if [ -s {input.vcf} ]; then "
            "  awk '{{print $2\"\\t\"$3\"\\t\"$6}}' {input.merged_nodups} | "
            "  sort -k1,1 -k2,2n > {params.out_dir}/fragments.txt && "
            "  extractHAIRS --hic 1 --VCF {input.vcf} "
            "  --fragments {params.out_dir}/fragments.txt "
            "  --out {params.out_dir}/fragment_file.txt 2>&1 | tee {log} && "
            "  HAPCUT2 --VCF {input.vcf} "
            "  --fragments {params.out_dir}/fragment_file.txt "
            "  --output {output.phased_vcf} "
            "  --hic 1 2>&1 | tee -a {log} && "
            "  python scripts/generate_phased_haplotypes.py "
            "  --vcf {output.phased_vcf} "
            "  --reference {input.mag} "
            "  --output1 {output.haplotype1} "
            "  --output2 {output.haplotype2} 2>&1 | tee -a {log}; "
            "else "
            "  echo 'No SNPs found, skipping phasing' | tee {log}; "
            "  touch {output.haplotype1} {output.haplotype2} {output.phased_vcf}; "
            "fi"

# ======================================================================
# 8. Optional: 3D-DNA scaffolding
# ======================================================================

if config.get("run_3d_dna", False):
    rule run_3d_dna:
        input:
            merged_nodups = "juicer/{sample}/{mag_name}/aligned/merged_nodups.txt",
            mag = "juicer/{sample}/references/{mag_name}/{mag_name}.fa",
            juicer_done = "juicer/{sample}/{mag_name}/juicer.done"
        output:
            assembly = "juicer/{sample}/{mag_name}/3d-dna/contigs.fasta",
            done = touch("juicer/{sample}/{mag_name}/3d-dna.done")
        benchmark: "logs/juicer/3d_dna/{sample}/{mag_name}.benchmark.txt"
        log: "logs/juicer/3d_dna/{sample}/{mag_name}.log"
        threads: config["n_threads"]
        params:
            juicer_path = config.get("juicer_path", "auto")
        conda: "envs/juicer.yaml"
        shell:
            "if [ \"{params.juicer_path}\" = \"auto\" ]; then "
            "  export JUICER_PATH=$(dirname $(dirname $(which juicer.sh))); "
            "else "
            "  export JUICER_PATH={params.juicer_path}; "
            "fi && "
            "OUT_DIR=juicer/{wildcards.sample}/{wildcards.mag_name}/3d-dna && "
            "mkdir -p $OUT_DIR && "
            "cd $OUT_DIR && "
            "bash $JUICER_PATH/3d-dna/run-asm-pipeline.sh "
            "-i 5000 --polisher-input-size 5000 "
            "$(realpath ../../references/{wildcards.mag_name}/{wildcards.mag_name}.fa) "
            "$(realpath ../../aligned/merged_nodups.txt) "
            "2>&1 | tee {log} && "
            "touch ../../3d-dna.done"

# ======================================================================
# 9. Completion marker
# ======================================================================

rule hic_maps_complete:
    input:
        mag_list = "juicer/{sample}/mag_list.txt",
        mags_done = "juicer/{sample}/mags_collected.done"
    output:
        done = touch("juicer/{sample}/all_hic_maps.done")
    benchmark: "logs/juicer/complete/{sample}.benchmark.txt"
    log: "logs/juicer/complete/{sample}.log"
    run:
        with open(input.mag_list) as f:
            mags = [line.strip() for line in f if line.strip()]

        missing = []
        for mag in mags:
            expected = f"juicer/{wildcards.sample}/{mag}/juicer.done"
            if not os.path.exists(expected):
                missing.append(mag)

        with open(output.done, 'w') as out:
            if missing:
                out.write(f"WARNING: {len(missing)}/{len(mags)} MAGs not processed: {', '.join(missing)}\n")
            else:
                out.write(f"All {len(mags)} MAGs processed successfully\n")
