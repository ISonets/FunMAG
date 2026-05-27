# Snakefile – Integrated WGS/Hi‑C fungal MAG reconstruction
# Final version with Hi‑C decontamination (human removal)
# Multi-binner → DASTool → BUSCO → Taxonomy → dRep → MetaWRAP → MAGPurify

import os

configfile: "config.yaml"
SAMPLES = config["samples"]

# Use filtered contigs after length filtering
ASSEMBLY_CONTIGS_FILTERED = "assembly/{sample}/contigs_filtered.fasta"

# ======================================================================
# Final targets
# ======================================================================
rule all:
    input:
        "results/merged_metaphlan_HiC.txt",
        "results/merged_metaphlan_WGS.txt",
        "results/merged_metaphlan_HiC_s.txt",
        "results/merged_metaphlan_WGS_s.txt",
        expand("results/{sample}/binning_stats.txt", sample=SAMPLES),
        expand("results/{sample}/median_cov_bin3c.tsv", sample=SAMPLES),
        expand("results/{sample}/median_cov_metabat.tsv", sample=SAMPLES),
        expand("results/{sample}/median_cov_maxbin2.tsv", sample=SAMPLES),
        expand("results/{sample}/fungal_bins_final.done", sample=SAMPLES),
        expand("results/{sample}/coverm_fungal_abundance.tsv", sample=SAMPLES),
        expand("results/{sample}/virsorter2_fungal.done", sample=SAMPLES),
        expand("results/{sample}/quast.done", sample=SAMPLES),
        expand("classification/{sample}/kaiju_report.html", sample=SAMPLES),
        "results/multiqc.done",

# ======================================================================
# 1. Human read removal (kneaddata) – BOTH WGS and Hi‑C
# ======================================================================

rule kneaddata_human:
    """
    Remove human reads using kneaddata (bowtie2 against hg38).
    Works for both WGS and Hi‑C reads via {read_type} wildcard.
    Output: decontaminated paired reads.
    """
    input:
        r1 = config["data_dir"] + "/{sample}_{read_type}_R1.fastq.gz",
        r2 = config["data_dir"] + "/{sample}_{read_type}_R2.fastq.gz"
    output:
        r1_clean = "decontam/{sample}_{read_type}_R1_clean.fastq.gz",
        r2_clean = "decontam/{sample}_{read_type}_R2_clean.fastq.gz",
    benchmark: "logs/kneaddata/{sample}_{read_type}.benchmark.txt"
    log: "logs/kneaddata/{sample}_{read_type}.log"
    threads: config["n_threads"]
    conda: "envs/kneaddata.yaml"
    params:
        db = config.get("kneaddata_db", "/path/to/kneaddata_human_db")
    shell:
        "kneaddata --input {input.r1} --input {input.r2} "
        "--reference-db {params.db} --output decontam/{wildcards.sample}_{wildcards.read_type} "
        "--threads {threads} --trimmomatic-options \"{config['trimmomatic_pe_opts']}\" "
        "--output-prefix {wildcards.sample}_{wildcards.read_type} 2>&1 | tee {log} && "
        "cp decontam/{wildcards.sample}_{wildcards.read_type}/{wildcards.sample}_{wildcards.read_type}_R1_paired.fastq.gz {output.r1_clean} && "
        "cp decontam/{wildcards.sample}_{wildcards.read_type}/{wildcards.sample}_{wildcards.read_type}_R2_paired.fastq.gz {output.r2_clean} || "
        "touch {output.r1_clean} {output.r2_clean}"

# ======================================================================
# 2. Quality control (on decontaminated reads)
# ======================================================================

rule fastqc:
    input:
        r1_wgs = "decontam/{sample}_WGS_R1_clean.fastq.gz",
        r2_wgs = "decontam/{sample}_WGS_R2_clean.fastq.gz",
        r1_hic = "decontam/{sample}_HiC_R1_clean.fastq.gz",
        r2_hic = "decontam/{sample}_HiC_R2_clean.fastq.gz"
    output:
        zip_wgs_r1 = "qc/fastqc/{sample}_WGS_R1_fastqc.zip",
        zip_wgs_r2 = "qc/fastqc/{sample}_WGS_R2_fastqc.zip",
        zip_hic_r1 = "qc/fastqc/{sample}_HiC_R1_fastqc.zip",
        zip_hic_r2 = "qc/fastqc/{sample}_HiC_R2_fastqc.zip"
    benchmark: "logs/fastqc/{sample}.benchmark.txt"
    log: "logs/fastqc/{sample}.log"
    threads: 4
    conda: "envs/bioinf_v1.yaml"
    shell: "fastqc -o qc/fastqc -f fastq -t {threads} {input} > {log} 2>&1"

rule seqprep:
    input:
        r1 = "decontam/{sample}_WGS_R1_clean.fastq.gz",
        r2 = "decontam/{sample}_WGS_R2_clean.fastq.gz"
    output:
        unmerged1 = "seqprep/{sample}/unmerged1.fastq.gz",
        unmerged2 = "seqprep/{sample}/unmerged2.fastq.gz",
        merged    = "seqprep/{sample}/merged.fastq.gz"
    params:
        adapt_a = config["adapters"]["A"],
        adapt_b = config["adapters"]["B"]
    benchmark: "logs/seqprep/{sample}.benchmark.txt"
    log: "logs/seqprep/{sample}.log"
    conda: "envs/hicmag_py37.yaml"
    shell:
        "SeqPrep -f {input.r1} -r {input.r2} -1 {output.unmerged1} -2 {output.unmerged2} "
        "-s {output.merged} -A {params.adapt_a} -B {params.adapt_b} 2>&1 | tee {log}"

rule trimmomatic_pe:
    input:
        unmerged1 = "seqprep/{sample}/unmerged1.fastq.gz",
        unmerged2 = "seqprep/{sample}/unmerged2.fastq.gz"
    output:
        pe1 = "trimmed/{sample}_WGS_R1_paired.fq.gz",
        ue1 = "trimmed/{sample}_WGS_R1_unpaired.fq.gz",
        pe2 = "trimmed/{sample}_WGS_R2_paired.fq.gz",
        ue2 = "trimmed/{sample}_WGS_R2_unpaired.fq.gz"
    benchmark: "logs/trimmomatic_pe/{sample}.benchmark.txt"
    log: "logs/trimmomatic_pe/{sample}.log"
    threads: config["n_threads"]
    conda: "envs/hic_mag.yaml"
    params: opts = config["trimmomatic_pe_opts"]
    shell:
        "trimmomatic PE -threads {threads} -phred33 {input} {output} {params.opts} 2>&1 | tee {log}"

rule trimmomatic_se:
    input: "seqprep/{sample}/merged.fastq.gz"
    output: "trimmed/{sample}_merged.fastq.gz"
    benchmark: "logs/trimmomatic_se/{sample}.benchmark.txt"
    log: "logs/trimmomatic_se/{sample}.log"
    threads: config["n_threads"]
    conda: "envs/hic_mag.yaml"
    params: opts = config["trimmomatic_se_opts"]
    shell:
        "trimmomatic SE -threads {threads} -phred33 {input} {output} {params.opts} 2>&1 | tee {log}"

# ======================================================================
# 3. Assembly (SPAdes meta or MEGAHIT)
# ======================================================================

if not config.get("use_megahit", False):
    rule spades:
        input:
            pe1 = "trimmed/{sample}_WGS_R1_paired.fq.gz",
            pe2 = "trimmed/{sample}_WGS_R2_paired.fq.gz",
            merged = "trimmed/{sample}_merged.fastq.gz"
        output: "assembly/{sample}/contigs.fasta"
        benchmark: "logs/spades/{sample}.benchmark.txt"
        log: "logs/spades/{sample}.log"
        threads: config["n_threads"]
        conda: "envs/hicmag_py37.yaml"
        params:
            kmers = config["spades_kmers"],
            mem = config["spades_mem"]
        shell:
            "spades.py --meta -1 {input.pe1} -2 {input.pe2} --merged {input.merged} "
            "-o assembly/{wildcards.sample} -t {threads} -k {params.kmers} -m {params.mem} 2>&1 | tee {log}"
else:
    rule megahit:
        input:
            pe1 = "trimmed/{sample}_WGS_R1_paired.fq.gz",
            pe2 = "trimmed/{sample}_WGS_R2_paired.fq.gz",
            merged = "trimmed/{sample}_merged.fastq.gz"
        output: "assembly/{sample}/megahit/final.contigs.fa"
        benchmark: "logs/megahit/{sample}.benchmark.txt"
        log: "logs/megahit/{sample}.log"
        threads: config["n_threads"]
        conda: "envs/hic_mag.yaml"
        shell:
            "megahit -t {threads} -1 {input.pe1} -2 {input.pe2} -r {input.merged} "
            "-o assembly/{wildcards.sample}/megahit 2>&1 | tee {log}"

# ======================================================================
# 4. Contig length filtering (≥1000 bp)
# ======================================================================

rule filter_contigs:
    input:
        fasta = "assembly/{sample}/megahit/final.contigs.fa" if config.get("use_megahit", False) 
               else "assembly/{sample}/contigs.fasta"
    output: ASSEMBLY_CONTIGS_FILTERED
    benchmark: "logs/filter_contigs/{sample}.benchmark.txt"
    log: "logs/filter_contigs/{sample}.log"
    conda: "envs/bioinf_v1.yaml"
    params: min_len = config["min_contig_length"]
    shell:
        "python scripts/filter_contigs_by_length.py {input.fasta} {output} --min-length {params.min_len} 2>&1 | tee {log}"

# ======================================================================
# 5. Index assembly & WGS alignment (PE only)
# ======================================================================

rule bwa_index:
    input: ASSEMBLY_CONTIGS_FILTERED
    output: touch(ASSEMBLY_CONTIGS_FILTERED + ".bwt")
    benchmark: "logs/bwa_index/{sample}.benchmark.txt"
    log: "logs/bwa_index/{sample}.log"
    conda: "envs/hic_mag.yaml"
    shell: "bwa index {input} 2>&1 | tee {log}"

rule bwa_mem_wgs:
    input:
        contigs = ASSEMBLY_CONTIGS_FILTERED,
        pe1 = "trimmed/{sample}_WGS_R1_paired.fq.gz",
        pe2 = "trimmed/{sample}_WGS_R2_paired.fq.gz"
    output: "binning/{sample}_wgs.sam"
    benchmark: "logs/bwa_mem_wgs/{sample}.benchmark.txt"
    log: "logs/bwa_mem_wgs/{sample}.log"
    threads: config["n_threads"]
    conda: "envs/hic_mag.yaml"
    shell: "bwa mem -t {threads} {input.contigs} {input.pe1} {input.pe2} > {output} 2>&1 | tee {log}"

rule samtools_sort_wgs:
    input: "binning/{sample}_wgs.sam"
    output: "binning/{sample}_wgs_sorted.bam"
    benchmark: "logs/samtools_sort_wgs/{sample}.benchmark.txt"
    log: "logs/samtools_sort_wgs/{sample}.log"
    threads: config["n_threads"]
    conda: "envs/bioinf_v1.yaml"
    shell: "samtools sort {input} -@ {threads} -o {output} 2>&1 | tee {log}"

# ======================================================================
# 6. Hi‑C processing: qc3C → bbduk → alignment → HiCzin → bin3C
# ======================================================================

if config.get("run_qc3c", True):
    rule qc3c:
        input:
            r1 = "decontam/{sample}_HiC_R1_clean.fastq.gz",
            r2 = "decontam/{sample}_HiC_R2_clean.fastq.gz"
        output:
            r1_qc = "qc3c/{sample}_HiC_R1.filt.fastq.gz",
            r2_qc = "qc3c/{sample}_HiC_R2.filt.fastq.gz"
        benchmark: "logs/qc3c/{sample}.benchmark.txt"
        log: "logs/qc3c/{sample}.log"
        conda: "envs/qc3c.yaml"
        params: extra = ""
        shell:
            "qc3C.py -1 {input.r1} -2 {input.r2} -o qc3c/{wildcards.sample} {params.extra} && "
            "mv qc3c/{wildcards.sample}/*_1.fastq.gz {output.r1_qc} && "
            "mv qc3c/{wildcards.sample}/*_2.fastq.gz {output.r2_qc} 2>&1 | tee {log}"

rule bbduk_hic:
    input:
        r1 = "qc3c/{sample}_HiC_R1.filt.fastq.gz" if config.get("run_qc3c", True) else "decontam/{sample}_HiC_R1_clean.fastq.gz",
        r2 = "qc3c/{sample}_HiC_R2.filt.fastq.gz" if config.get("run_qc3c", True) else "decontam/{sample}_HiC_R2_clean.fastq.gz"
    output: "hic_trimmed/{sample}_paired.fastq.gz"
    benchmark: "logs/bbduk_hic/{sample}.benchmark.txt"
    log: "logs/bbduk_hic/{sample}.log"
    conda: "envs/hic_mag.yaml"
    shell:
        "bbduk.sh in1={input.r1} in2={input.r2} k=23 hdist=1 mink=11 "
        "ktrim=r tpe tbo ftm=5 qtrim=r trimq=10 out={output} 2>&1 | tee {log}"

rule bwa_mem_hic:
    input:
        contigs = ASSEMBLY_CONTIGS_FILTERED,
        hic = "hic_trimmed/{sample}_paired.fastq.gz"
    output: "binning/{sample}_hic_unsorted.bam"
    benchmark: "logs/bwa_mem_hic/{sample}.benchmark.txt"
    log: "logs/bwa_mem_hic/{sample}.log"
    threads: config["n_threads"]
    conda: "envs/hic_mag.yaml"
    shell:
        "bwa mem -5SP -t {threads} {input.contigs} {input.hic} | "
        "samtools view -F 0x904 -bS -o {output} - 2>&1 | tee {log}"

rule samtools_sort_hic_name:
    input: "binning/{sample}_hic_unsorted.bam"
    output: "binning/{sample}_hic_sorted_n.bam"
    benchmark: "logs/samtools_sort_hic_name/{sample}.benchmark.txt"
    log: "logs/samtools_sort_hic_name/{sample}.log"
    threads: config["n_threads"]
    conda: "envs/bioinf_v1.yaml"
    shell: "samtools sort -n -@ {threads} -o {output} {input} 2>&1 | tee {log}"

rule bin3c_mkmap:
    input:
        contigs = ASSEMBLY_CONTIGS_FILTERED,
        bam = "binning/{sample}_hic_sorted_n.bam"
    output: "bin3c/{sample}/contact_map.p.gz"
    benchmark: "logs/bin3c_mkmap/{sample}.benchmark.txt"
    log: "logs/bin3c_mkmap/{sample}.log"
    params:
        enzyme = config["restriction_enzyme"],
        resolution = config["bin3c_resolution"]
    conda: "envs/bin3c.yaml"
    shell:
        "bin3C.py mkmap -e {params.enzyme} -r {params.resolution} "
        "-v {input.contigs} {input.bam} bin3c/{wildcards.sample} 2>&1 | tee {log}"

rule hiczin_normalize:
    input: "bin3c/{sample}/contact_map.p.gz"
    output: "bin3c/{sample}/contact_map_hiczin.p.gz"
    benchmark: "logs/hiczin/{sample}.benchmark.txt"
    log: "logs/hiczin/{sample}.log"
    params:
        iterations = config["hiczin_iterations"],
        filter_threshold = config["hiczin_filter_threshold"]
    conda: "envs/hiczin.yaml"
    shell:
        "HiCzin.py -i {input} -o {output} "
        "--iterations {params.iterations} --filter-threshold {params.filter_threshold} 2>&1 | tee {log}"

rule bin3c_cluster:
    input: "bin3c/{sample}/contact_map_hiczin.p.gz"
    output: directory("bin3c/{sample}/fasta")
    benchmark: "logs/bin3c_cluster/{sample}.benchmark.txt"
    log: "logs/bin3c_cluster/{sample}.log"
    conda: "envs/bin3c.yaml"
    shell: "bin3C.py cluster --only-large -v {input} bin3c/{wildcards.sample} 2>&1 | tee {log}"

# ======================================================================
# 7. WGS binning: MetaBAT2 + MaxBin2
# ======================================================================

rule metabat2:
    input:
        contigs = ASSEMBLY_CONTIGS_FILTERED,
        bam = "binning/{sample}_wgs_sorted.bam"
    output: directory("metabat/{sample}_bins")
    benchmark: "logs/metabat2/{sample}.benchmark.txt"
    log: "logs/metabat2/{sample}.log"
    threads: config["n_threads"]
    params: min_contig = config["metabat2_min_contig"]
    conda: "envs/bioinf_v1.yaml"
    shell:
        "metabat2 -i {input.contigs} {input.bam} -o metabat/{wildcards.sample}_bins/bin "
        "-v -t {threads} --minContig {params.min_contig} 2>&1 | tee {log}"

rule jgi_summarize_depth:
    input: "binning/{sample}_wgs_sorted.bam"
    output: "metabat/{sample}_depth.tsv"
    benchmark: "logs/jgi_summarize_depth/{sample}.benchmark.txt"
    log: "logs/jgi_summarize_depth/{sample}.log"
    conda: "envs/bioinf_v1.yaml"
    shell: "jgi_summarize_bam_contig_depths --outputDepth {output} {input} 2>&1 | tee {log}"

rule maxbin2:
    input:
        contigs = ASSEMBLY_CONTIGS_FILTERED,
        depth = "metabat/{sample}_depth.tsv"
    output:
        bins = directory("maxbin2/{sample}_bins"),
        summary = "maxbin2/{sample}_maxbin2.summary"
    benchmark: "logs/maxbin2/{sample}.benchmark.txt"
    log: "logs/maxbin2/{sample}.log"
    threads: config["n_threads"]
    params:
        out_prefix = "maxbin2/{sample}_bins/bin",
        prob_threshold = config["maxbin2_prob_threshold"],
        min_contig = config["maxbin2_min_contig"]
    conda: "envs/maxbin2.yaml"
    shell:
        "cut -f1,4- {input.depth} | sed 's/contigName/Contig\\t{wildcards.sample}/' > maxbin2/{wildcards.sample}_abund.tsv && "
        "run_MaxBin.pl -contig {input.contigs} -abund maxbin2/{wildcards.sample}_abund.tsv "
        "-out {params.out_prefix} -thread {threads} -min_contig_length {params.min_contig} "
        "-max_iteration 50 -prob_threshold {params.prob_threshold} 2>&1 | tee {log} && "
        "mkdir -p {output.bins} && "
        "cp {params.out_prefix}.*.fasta {output.bins}/"

# ======================================================================
# 8. DASTool consensus binning (MetaBAT2 + MaxBin2 + bin3C)
# ======================================================================

rule generate_hic_contig_list:
    input: "bin3c/{sample}/fasta"
    output: "dastool/{sample}/hic_contig_list.txt"
    benchmark: "logs/generate_hic_contig_list/{sample}.benchmark.txt"
    log: "logs/generate_hic_contig_list/{sample}.log"
    run:
        with open(output[0], 'w') as out:
            for fna in sorted(os.listdir(input[0])):
                if not fna.endswith('.fna'): continue
                with open(os.path.join(input[0], fna)) as f:
                    contigs = [line[1:].strip().split()[0] for line in f if line.startswith('>')]
                if contigs:
                    out.write(' '.join(contigs) + '\n')

rule generate_wgs_contig_lists:
    input:
        metabat = "metabat/{sample}_bins",
        maxbin2 = "maxbin2/{sample}_bins"
    output:
        metabat_lst = "dastool/{sample}/metabat_contig_list.txt",
        maxbin2_lst = "dastool/{sample}/maxbin2_contig_list.txt"
    benchmark: "logs/generate_wgs_contig_lists/{sample}.benchmark.txt"
    log: "logs/generate_wgs_contig_lists/{sample}.log"
    run:
        for bin_dir, out_file in [(input.metabat, output.metabat_lst), (input.maxbin2, output.maxbin2_lst)]:
            with open(out_file, 'w') as out:
                if os.path.isdir(bin_dir):
                    for fa in sorted(os.listdir(bin_dir)):
                        if not fa.endswith(('.fa', '.fasta')): continue
                        with open(os.path.join(bin_dir, fa)) as f:
                            contigs = [line[1:].strip().split()[0] for line in f if line.startswith('>')]
                        if contigs:
                            out.write(' '.join(contigs) + '\n')

rule prep_tsv_dastool:
    input:
        metabat_lst = "dastool/{sample}/metabat_contig_list.txt",
        maxbin2_lst = "dastool/{sample}/maxbin2_contig_list.txt",
        hic_lst = "dastool/{sample}/hic_contig_list.txt"
    output:
        metabat_tsv = "dastool/{sample}/metabat_contig2bin.tsv",
        maxbin2_tsv = "dastool/{sample}/maxbin2_contig2bin.tsv",
        hic_tsv = "dastool/{sample}/hic_contig2bin.tsv"
    benchmark: "logs/prep_tsv_dastool/{sample}.benchmark.txt"
    log: "logs/prep_tsv_dastool/{sample}.log"
    conda: "envs/bioinf_v1.yaml"
    shell:
        "python scripts/prep_tsv_for_dastool_multi.py "
        "{input.metabat_lst} metabat {output.metabat_tsv} "
        "{input.maxbin2_lst} maxbin2 {output.maxbin2_tsv} "
        "{input.hic_lst} bin3c {output.hic_tsv} 2>&1 | tee {log}"

if config.get("run_dastool", True):
    rule dastool:
        input:
            contigs = ASSEMBLY_CONTIGS_FILTERED,
            metabat_tsv = "dastool/{sample}/metabat_contig2bin.tsv",
            maxbin2_tsv = "dastool/{sample}/maxbin2_contig2bin.tsv",
            hic_tsv = "dastool/{sample}/hic_contig2bin.tsv"
        output: directory("dastool/{sample}/DASTool_bins")
        benchmark: "logs/dastool/{sample}.benchmark.txt"
        log: "logs/dastool/{sample}.log"
        threads: config["n_threads"]
        params:
            out_dir = "dastool/{sample}",
            engine = config["dastool_search_engine"]
        conda: "envs/dastool.yaml"
        shell:
            "DAS_Tool -i {input.metabat_tsv},{input.maxbin2_tsv},{input.hic_tsv} "
            "-c {input.contigs} -o {params.out_dir}/DASTool --write_bins 1 "
            "--search_engine {params.engine} -t {threads} 2>&1 | tee {log}"

    rule krakenuniq_dastool:
        input: "dastool/{sample}/DASTool_bins"
        output: "classification/{sample}/krakenuniq_dastool.report"
        benchmark: "logs/krakenuniq_dastool/{sample}.benchmark.txt"
        log: "logs/krakenuniq_dastool/{sample}.log"
        threads: config["n_threads"]
        params: db = config["krakenuniq_db"]
        conda: "envs/krakenuniq.yaml"
        shell:
            "cat {input}/*.fa > tmp_dastool_all.fa; "
            "krakenuniq --db {params.db} --threads {threads} "
            "--report-file {output} --output /dev/null tmp_dastool_all.fa 2>&1 | tee {log}"

# ======================================================================
# 9. Collect ALL bins from all binners + DASTool
# ======================================================================

rule collect_all_bins:
    input:
        metabat = "metabat/{sample}_bins",
        maxbin2 = "maxbin2/{sample}_bins",
        bin3c = "bin3c/{sample}/fasta",
        dastool = "dastool/{sample}/DASTool_bins" if config.get("run_dastool", True) else []
    output: "classification/{sample}/all_bins.lst"
    benchmark: "logs/collect_all_bins/{sample}.benchmark.txt"
    log: "logs/collect_all_bins/{sample}.log"
    run:
        bins = []
        for d in [input.metabat, input.maxbin2, input.bin3c]:
            if os.path.isdir(d):
                for f in sorted(os.listdir(d)):
                    if f.endswith(('.fna', '.fa', '.fasta')):
                        bins.append(os.path.abspath(os.path.join(d, f)))
        dastool_dir = f"dastool/{wildcards.sample}/DASTool_bins"
        if os.path.isdir(dastool_dir):
            for f in sorted(os.listdir(dastool_dir)):
                if f.endswith(('.fna', '.fa', '.fasta')):
                    bins.append(os.path.abspath(os.path.join(dastool_dir, f)))
        with open(output[0], 'w') as out:
            for b in bins:
                out.write(b + "\n")

# ======================================================================
# 10. BUSCO on ALL bins (quality filtering)
# ======================================================================

rule busco_all_bins:
    input:
        lst = "classification/{sample}/all_bins.lst"
    output:
        tsv = "busco/{sample}/busco_results.tsv",
        filtered_lst = "classification/{sample}/bins_busco_pass.lst"
    benchmark: "logs/busco/{sample}.benchmark.txt"
    log: "logs/busco/{sample}.log"
    threads: config["n_threads"]
    params:
        lineage = config.get("busco_lineage", "auto"),
        min_completeness = config["busco_min_completeness"],
        max_contamination = config["busco_max_contamination"],
        mode = config.get("busco_mode", "comprehensive")
    conda: "envs/busco.yaml"
    shell:
        "python scripts/run_busco_filter.py --bin-list {input.lst} "
        "--output {output.tsv} --filtered {output.filtered_lst} "
        "--lineage {params.lineage} --threads {threads} "
        "--min-completeness {params.min_completeness} "
        "--max-contamination {params.max_contamination} "
        "--mode {params.mode} 2>&1 | tee {log}"

# ======================================================================
# 11. Community taxonomic profiling (MetaPhlAn, MiCoP, Kraken2)
# ======================================================================

rule metaphlan_wgs:
    input: "trimmed/{sample}_merged.fastq.gz"
    output:
        profile = "profiling/{sample}/metaphlan_wgs.txt",
        bz2 = "profiling/{sample}/metaphlan_wgs.bowtie2.bz2"
    benchmark: "logs/metaphlan_wgs/{sample}.benchmark.txt"
    log: "logs/metaphlan_wgs/{sample}.log"
    threads: config["n_threads"]
    conda: "envs/hic_mag.yaml"
    shell:
        "metaphlan {input} --bowtie2out {output.bz2} --nproc {threads} "
        "--input_type fastq -o {output.profile} 2>&1 | tee {log}"

rule metaphlan_wgs_species:
    input: "profiling/{sample}/metaphlan_wgs.bowtie2.bz2"
    output: "profiling/{sample}/metaphlan_wgs_species.txt"
    benchmark: "logs/metaphlan_wgs_species/{sample}.benchmark.txt"
    log: "logs/metaphlan_wgs_species/{sample}.log"
    threads: config["n_threads"]
    conda: "envs/hic_mag.yaml"
    shell:
        "metaphlan {input} --nproc {threads} --input_type bowtie2out "
        "-o {output} --tax_lev 's' 2>&1 | tee {log}"

rule metaphlan_hic:
    input: "hic_trimmed/{sample}_paired.fastq.gz"
    output:
        profile = "profiling/{sample}/metaphlan_hic.txt",
        bz2 = "profiling/{sample}/metaphlan_hic.bowtie2.bz2"
    benchmark: "logs/metaphlan_hic/{sample}.benchmark.txt"
    log: "logs/metaphlan_hic/{sample}.log"
    threads: config["n_threads"]
    conda: "envs/hic_mag.yaml"
    shell:
        "metaphlan {input} --bowtie2out {output.bz2} --nproc {threads} "
        "--input_type fastq -o {output.profile} 2>&1 | tee {log}"

rule metaphlan_hic_species:
    input: "profiling/{sample}/metaphlan_hic.bowtie2.bz2"
    output: "profiling/{sample}/metaphlan_hic_species.txt"
    benchmark: "logs/metaphlan_hic_species/{sample}.benchmark.txt"
    log: "logs/metaphlan_hic_species/{sample}.log"
    threads: config["n_threads"]
    conda: "envs/hic_mag.yaml"
    shell:
        "metaphlan {input} --nproc {threads} --input_type bowtie2out "
        "-o {output} --tax_lev 's' 2>&1 | tee {log}"

rule micop_fungi_community:
    input: "trimmed/{sample}_merged.fastq.gz"
    output:
        sam = "profiling/{sample}/micop_fungi.sam",
        abund = "profiling/{sample}/micop_fungi.txt"
    benchmark: "logs/micop_fungi/{sample}.benchmark.txt"
    log: "logs/micop_fungi/{sample}.log"
    conda: "envs/bioinf_v1.yaml"
    shell:
        "python tools/MiCoP/run-bwa.py {input} --fungi --output {output.sam} && "
        "python tools/MiCoP/compute-abundances.py {output.sam} --fungi --output {output.abund} 2>&1 | tee {log}"

rule kraken2_reads:
    input:
        r1 = "decontam/{sample}_WGS_R1_clean.fastq.gz",
        r2 = "decontam/{sample}_WGS_R2_clean.fastq.gz"
    output: "profiling/{sample}/kraken2_report.txt"
    benchmark: "logs/kraken2_reads/{sample}.benchmark.txt"
    log: "logs/kraken2_reads/{sample}.log"
    threads: config["n_threads"]
    params: db = config["kraken2_db"]
    conda: "envs/hic_mag.yaml"
    shell:
        "kraken2 --db {params.db} --threads {threads} --report {output} "
        "--paired {input.r1} {input.r2} > /dev/null 2>&1 | tee {log}"

# ======================================================================
# 12. Bin taxonomy: BLAST nr, Kaiju, MetaPhlAn, MiCoP (on BUSCO-passed bins)
# ======================================================================

rule collect_busco_bins:
    input: "classification/{sample}/bins_busco_pass.lst"
    output: "classification/{sample}/busco_bins_for_taxonomy.lst"
    shell: "cp {input} {output}"

rule blast_nr:
    input: "classification/{sample}/busco_bins_for_taxonomy.lst"
    output: "classification/{sample}/blast_nr_top10.csv"
    benchmark: "logs/blast_nr/{sample}.benchmark.txt"
    log: "logs/blast_nr/{sample}.log"
    threads: config["n_threads"]
    params:
        db = config["blast_nr_db"],
        outfmt = "6 qseqid sseqid pident length stitle"
    conda: "envs/blast.yaml"
    shell:
        "while read bin; do "
        "  blastn -db {params.db} -query $bin -max_target_seqs 10 "
        "    -outfmt \"{params.outfmt}\" -num_threads {threads} -out tmp.tsv; "
        "  echo -n \"$bin\t\"; paste -s tmp.tsv; "
        "done < {input} > {output} 2>&1 | tee {log}"

rule kaiju_bins:
    input: "classification/{sample}/busco_bins_for_taxonomy.lst"
    output: "classification/{sample}/kaiju.out"
    benchmark: "logs/kaiju_bins/{sample}.benchmark.txt"
    log: "logs/kaiju_bins/{sample}.log"
    threads: config["n_threads"]
    params:
        db = config["kaiju_db"],
        nodes = config["kaiju_nodes"]
    conda: "envs/kaiju.yaml"
    shell:
        "while read bin; do "
        "  kaiju -t {params.nodes} -f {params.db} -i $bin -o tmp.kaiju -z {threads} -v; "
        "  kaiju2table -t {params.nodes} -n {params.db} -r genus -o tmp.kaiju.table tmp.kaiju; "
        "  echo -n \"$bin\t\"; awk 'NR>1' tmp.kaiju.table | tr '\n' ';'; echo; "
        "done < {input} > {output} 2>&1 | tee {log}"

rule metaphlan_bins:
    input: "classification/{sample}/busco_bins_for_taxonomy.lst"
    output: "classification/{sample}/metaphlan_bins.csv"
    benchmark: "logs/metaphlan_bins/{sample}.benchmark.txt"
    log: "logs/metaphlan_bins/{sample}.log"
    threads: config["n_threads"]
    conda: "envs/hic_mag.yaml"
    shell:
        "echo 'bin,taxonomy' > {output}; "
        "while read bin; do "
        "  metaphlan $bin --input_type fasta --nproc {threads} -o tmp_mpa.txt; "
        "  echo -n \"$bin,\"; "
        "  grep -E '(k__|p__|c__)' tmp_mpa.txt | paste -s -d ';'; "
        "done < {input} >> {output} 2>&1 | tee {log}"

rule micop_bins_fungi:
    input: "classification/{sample}/busco_bins_for_taxonomy.lst"
    output: "classification/{sample}/micop_bins_fungi.csv"
    benchmark: "logs/micop_bins_fungi/{sample}.benchmark.txt"
    log: "logs/micop_bins_fungi/{sample}.log"
    conda: "envs/bioinf_v1.yaml"
    shell:
        "echo 'bin,fungal_reads,total_reads,fraction' > {output}; "
        "while read bin; do "
        "  python tools/MiCoP/run-bwa.py $bin --fungi --output tmp_fungi.sam; "
        "  python tools/MiCoP/compute-abundances.py tmp_fungi.sam --fungi --output tmp_fungi.txt; "
        "  reads=$(grep -w 'fungi' tmp_fungi.txt | cut -f2); total=$(grep -c '^>' $bin); "
        "  [[ -z $reads ]] && reads=0; "
        "  frac=$(echo \"scale=4; $reads/$total\" | bc); "
        "  echo \"$bin,$reads,$total,$frac\" >> {output}; "
        "done < {input} 2>&1 | tee {log}"

# ======================================================================
# 13. Fungal bin identification (consensus + BUSCO)
# ======================================================================

rule classify_and_filter_fungal:
    input:
        blast = "classification/{sample}/blast_nr_top10.csv",
        kaiju = "classification/{sample}/kaiju.out",
        metaphlan = "classification/{sample}/metaphlan_bins.csv",
        micop = "classification/{sample}/micop_bins_fungi.csv",
        busco_filtered = "classification/{sample}/bins_busco_pass.lst",
        dastool_kraken = "classification/{sample}/krakenuniq_dastool.report" if config.get("run_dastool", False) else [],
    output:
        done = touch("results/{sample}/fungal_bins_identified.done"),
        lst = "results/{sample}/fungal_bins_identified.lst",
        dir = directory("results/{sample}/fungal_bins_identified")
    benchmark: "logs/classify_and_filter/{sample}.benchmark.txt"
    log: "logs/classify_and_filter/{sample}.log"
    params:
        min_frac = config["min_fungal_fraction"],
        out_dir = "results/{sample}/fungal_bins_identified"
    conda: "envs/bioinf_v1.yaml"
    script: "scripts/classify_and_filter_bins.py"

# ======================================================================
# 14. dRep dereplication (deduplication of fungal bins)
# ======================================================================

rule drep_dereplicate:
    input:
        fungal_lst = "results/{sample}/fungal_bins_identified.lst",
        busco_tsv = "busco/{sample}/busco_results.tsv"
    output:
        done = touch("results/{sample}/drep_dereplicate.done"),
        representatives = "results/{sample}/fungal_bins_drep/fungal_representatives.lst",
        dir = directory("results/{sample}/fungal_bins_drep")
    benchmark: "logs/drep/{sample}.benchmark.txt"
    log: "logs/drep/{sample}.log"
    threads: config["n_threads"]
    params:
        out_dir = "results/{sample}/fungal_bins_drep",
        ani_threshold = config["drep_ani_threshold"],
        coverage = config["drep_coverage"]
    conda: "envs/drep.yaml"
    shell:
        "mkdir -p {params.out_dir} && "
        "python scripts/run_drep.py --bin-list {input.fungal_lst} "
        "--busco {input.busco_tsv} --output {params.out_dir} "
        "--ani {params.ani_threshold} --coverage {params.coverage} "
        "--threads {threads} 2>&1 | tee {log} && "
        "touch {output.done}"

# ======================================================================
# 15. MetaWRAP reassembly (only on dRep representatives)
# ======================================================================

rule metawrap_reassemble:
    input:
        reps_lst = "results/{sample}/fungal_bins_drep/fungal_representatives.lst",
        r1 = "trimmed/{sample}_WGS_R1_paired.fq.gz",
        r2 = "trimmed/{sample}_WGS_R2_paired.fq.gz"
    output:
        done = touch("results/{sample}/metawrap_reassemble.done"),
        dir = directory("results/{sample}/fungal_bins_reassembled")
    benchmark: "logs/metawrap_reassemble/{sample}.benchmark.txt"
    log: "logs/metawrap_reassemble/{sample}.log"
    threads: config["n_threads"]
    params:
        out_dir = "results/{sample}/fungal_bins_reassembled"
    conda: "envs/metawrap.yaml"
    shell:
        "mkdir -p {params.out_dir}/input_bins && "
        "python scripts/copy_bins_from_list.py {input.reps_lst} {params.out_dir}/input_bins && "
        "metaWRAP reassemble_bins -o {params.out_dir}/reassembly "
        "-b {params.out_dir}/input_bins -1 {input.r1} -2 {input.r2} "
        "-t {threads} --strict-cut-off 2 --permissive-cut-off 5 "
        "-c {config[busco_min_completeness]} -x {config[busco_max_contamination]} 2>&1 | tee {log} && "
        "touch {output.done}"

# ======================================================================
# 16. Re-run BUSCO on reassembled bins (verify improvement)
# ======================================================================

rule busco_reassembled:
    input: "results/{sample}/fungal_bins_reassembled/reassembly/reassembled_bins"
    output:
        tsv = "results/{sample}/busco_reassembled.tsv",
        filtered = "results/{sample}/fungal_bins_final.lst"
    benchmark: "logs/busco_reassembled/{sample}.benchmark.txt"
    log: "logs/busco_reassembled/{sample}.log"
    threads: config["n_threads"]
    params:
        lineage = config["busco_lineage"],
        min_completeness = config["busco_min_completeness"],
        max_contamination = config["busco_max_contamination"]
    conda: "envs/busco.yaml"
    shell:
        "ls {input}/*.fa 2>/dev/null > reassembled_bins.lst && "
        "if [ -s reassembled_bins.lst ]; then "
        "  python scripts/run_busco_filter.py --bin-list reassembled_bins.lst "
        "  --output {output.tsv} --filtered {output.filtered} "
        "  --lineage {params.lineage} --threads {threads} "
        "  --min-completeness {params.min_completeness} "
        "  --max-contamination {params.max_contamination} 2>&1 | tee {log}; "
        "else "
        "  echo 'No reassembled bins found' | tee {log}; "
        "  touch {output.tsv} {output.filtered}; "
        "fi"

# ======================================================================
# 17. MAGPurify (only on final BUSCO-passed reassembled bins)
# ======================================================================

rule magpurify:
    input:
        lst = "results/{sample}/fungal_bins_final.lst",
        bam = "binning/{sample}_wgs_sorted.bam"
    output:
        done = touch("results/{sample}/magpurify.done"),
        dir = directory("results/{sample}/fungal_bins_magpurify")
    benchmark: "logs/magpurify/{sample}.benchmark.txt"
    log: "logs/magpurify/{sample}.log"
    threads: config["n_threads"]
    params:
        out_dir = "results/{sample}/fungal_bins_magpurify"
    conda: "envs/magpurify.yaml"
    shell:
        "mkdir -p {params.out_dir} && "
        "if [ -s {input.lst} ]; then "
        "  while read bin; do "
        "    name=$(basename $bin .fa); "
        "    name=$(basename $name .fna); "
        "    name=$(basename $name .fasta); "
        "    magpurify tetra-freq $bin {params.out_dir} 2>/dev/null; "
        "    magpurify gc-content $bin {params.out_dir} 2>/dev/null; "
        "    if [ -f {input.bam} ]; then "
        "      magpurify coverage $bin {params.out_dir} --bam {input.bam} 2>/dev/null; "
        "    fi; "
        "    magpurify clean-bin $bin {params.out_dir}/${{name}}_clean.fa 2>/dev/null || "
        "    cp $bin {params.out_dir}/${{name}}_clean.fa; "
        "  done < {input.lst}; "
        "fi && "
        "touch {output.done}"

# ======================================================================
# 18. Final fungal bins done marker
# ======================================================================

rule fungal_bins_final:
    input:
        magpurify = "results/{sample}/magpurify.done",
        busco_final = "results/{sample}/busco_reassembled.tsv"
    output: touch("results/{sample}/fungal_bins_final.done")
    shell: "touch {output}"

# ======================================================================
# 19. QUAST reports for final fungal bins
# ======================================================================

rule quast_fungal_bins:
    input:
        lst = "results/{sample}/fungal_bins_final.lst",
        magpurify_done = "results/{sample}/magpurify.done"
    output:
        html = "results/{sample}/quast/report.html",
        done = touch("results/{sample}/quast.done")
    benchmark: "logs/quast/{sample}.benchmark.txt"
    log: "logs/quast/{sample}.log"
    threads: config["n_threads"]
    params:
        out_dir = "results/{sample}/quast"
    conda: "envs/quast.yaml"
    shell:
        "mkdir -p {params.out_dir} && "
        "if [ -s {input.lst} ]; then "
        "  quast.py $(cat {input.lst}) "
        "  -o {params.out_dir} "
        "  --threads {threads} "
        "  --fungus "
        "  --gene-finding "
        "  --rna-finding "
        "  --conserved-genes-finding "
        "  --circos "
        "  --html-report {output.html} "
        "  2>&1 | tee {log}; "
        "else "
        "  echo 'No fungal bins found for QUAST' | tee {log}; "
        "  touch {output.html}; "
        "fi && "
        "touch {output.done}"

# ======================================================================
# 20. Kaiju HTML report
# ======================================================================

rule kaiju_html_report:
    input:
        kaiju_out = "classification/{sample}/kaiju.out",
        busco_bins = "classification/{sample}/busco_bins_for_taxonomy.lst"
    output:
        html = "classification/{sample}/kaiju_report.html"
    benchmark: "logs/kaiju_html/{sample}.benchmark.txt"
    log: "logs/kaiju_html/{sample}.log"
    threads: 4
    conda: "envs/kaiju.yaml"
    params:
        nodes = config["kaiju_nodes"],
        names = config.get("kaiju_names", config["kaiju_nodes"].replace("nodes.dmp", "names.dmp"))
    shell:
        "if [ -s {input.kaiju_out} ]; then "
        "  kaiju2krona -t {params.nodes} -n {params.names} "
        "  -i {input.kaiju_out} -o kaiju_krona.txt 2>&1 | tee {log} && "
        "  ktImportText -o {output.html} kaiju_krona.txt 2>&1 | tee -a {log}; "
        "else "
        "  echo 'No Kaiju output to visualize' | tee {log}; "
        "  touch {output.html}; "
        "fi"

# ======================================================================
# 21. VirSorter2 (only on final fungal MAGs)
# ======================================================================

rule virsorter2:
    input:
        lst = "results/{sample}/fungal_bins_final.lst"
    output: touch("results/{sample}/virsorter2_fungal.done")
    benchmark: "logs/virsorter2/{sample}.benchmark.txt"
    log: "logs/virsorter2/{sample}.log"
    threads: config["n_threads"]
    params:
        out_dir = "results/{sample}/virsorter2",
        db = config.get("virsorter2_db")
    conda: "envs/virsorter2.yaml"
    shell:
        "mkdir -p {params.out_dir} && "
        "while read bin; do "
        "  name=$(basename $bin .fa); "
        "  virsorter run --seqfile $bin --include-groups dsDNAphage,RNA,ssDNA "
        "  --min-length 5000 --min-score 0.5 -j {threads} -w {params.out_dir}/$name; "
        "done < {input.lst} 2>&1 | tee {log} && touch {output}"

# ======================================================================
# 22. CoverM abundance (only on final fungal MAGs)
# ======================================================================

rule coverm_abundance:
    input:
        lst = "results/{sample}/fungal_bins_final.lst",
        r1 = "trimmed/{sample}_WGS_R1_paired.fq.gz",
        r2 = "trimmed/{sample}_WGS_R2_paired.fq.gz"
    output: "results/{sample}/coverm_fungal_abundance.tsv"
    benchmark: "logs/coverm/{sample}.benchmark.txt"
    log: "logs/coverm/{sample}.log"
    threads: config["n_threads"]
    params:
        out_dir = "results/{sample}/coverm",
        min_frac = config["coverm_min_covered_fraction"]
    conda: "envs/coverm.yaml"
    shell:
        "mkdir -p {params.out_dir} && "
        "if [ -s {input.lst} ]; then "
        "  cat $(cat {input.lst}) > {params.out_dir}/all_fungal_mags.fa && "
        "  coverm genome -d {params.out_dir}/all_fungal_mags.fa "
        "  -x fa -1 {input.r1} -2 {input.r2} "
        "  -t {threads} -m relative_abundance --min-covered-fraction {params.min_frac} "
        "  --output-file {output} 2>&1 | tee {log}; "
        "else "
        "  echo -e 'genome\\trelative_abundance' > {output}; "
        "fi"

# ======================================================================
# 23. Median coverage (Python scripts)
# ======================================================================

rule median_cov_bin3c:
    input:
        depth = "metabat/{sample}_depth.tsv",
        bins = "bin3c/{sample}/fasta"
    output: "results/{sample}/median_cov_bin3c.tsv"
    benchmark: "logs/median_cov_bin3c/{sample}.benchmark.txt"
    log: "logs/median_cov_bin3c/{sample}.log"
    conda: "envs/bioinf_v1.yaml"
    shell: "python scripts/median_cov_bin3c.py {input.depth} {input.bins} {output} 2>&1 | tee {log}"

rule median_cov_metabat:
    input:
        depth = "metabat/{sample}_depth.tsv",
        bins = "metabat/{sample}_bins"
    output: "results/{sample}/median_cov_metabat.tsv"
    benchmark: "logs/median_cov_metabat/{sample}.benchmark.txt"
    log: "logs/median_cov_metabat/{sample}.log"
    conda: "envs/bioinf_v1.yaml"
    shell: "python scripts/median_cov_metabat.py {input.depth} {input.bins} {output} 2>&1 | tee {log}"

rule median_cov_maxbin2:
    input:
        depth = "metabat/{sample}_depth.tsv",
        bins = "maxbin2/{sample}_bins"
    output: "results/{sample}/median_cov_maxbin2.tsv"
    benchmark: "logs/median_cov_maxbin2/{sample}.benchmark.txt"
    log: "logs/median_cov_maxbin2/{sample}.log"
    conda: "envs/bioinf_v1.yaml"
    shell: "python scripts/median_cov_maxbin2.py {input.depth} {input.bins} {output} 2>&1 | tee {log}"

# ======================================================================
# 24. Merge MetaPhlAn tables
# ======================================================================

rule merge_metaphlan_hic:
    input: expand("profiling/{sample}/metaphlan_hic.txt", sample=SAMPLES)
    output: "results/merged_metaphlan_HiC.txt"
    benchmark: "logs/merge_metaphlan_hic/benchmark.txt"
    log: "logs/merge_metaphlan_hic/log.txt"
    conda: "envs/hic_mag.yaml"
    shell: "merge_metaphlan_tables.py {input} > {output} 2>&1 | tee {log}"

rule merge_metaphlan_wgs:
    input: expand("profiling/{sample}/metaphlan_wgs.txt", sample=SAMPLES)
    output: "results/merged_metaphlan_WGS.txt"
    benchmark: "logs/merge_metaphlan_wgs/benchmark.txt"
    log: "logs/merge_metaphlan_wgs/log.txt"
    conda: "envs/hic_mag.yaml"
    shell: "merge_metaphlan_tables.py {input} > {output} 2>&1 | tee {log}"

rule merge_metaphlan_hic_s:
    input: expand("profiling/{sample}/metaphlan_hic_species.txt", sample=SAMPLES)
    output: "results/merged_metaphlan_HiC_s.txt"
    benchmark: "logs/merge_metaphlan_hic_s/benchmark.txt"
    log: "logs/merge_metaphlan_hic_s/log.txt"
    conda: "envs/hic_mag.yaml"
    shell: "merge_metaphlan_tables.py {input} > {output} 2>&1 | tee {log}"

rule merge_metaphlan_wgs_s:
    input: expand("profiling/{sample}/metaphlan_wgs_species.txt", sample=SAMPLES)
    output: "results/merged_metaphlan_WGS_s.txt"
    benchmark: "logs/merge_metaphlan_wgs_s/benchmark.txt"
    log: "logs/merge_metaphlan_wgs_s/log.txt"
    conda: "envs/hic_mag.yaml"
    shell: "merge_metaphlan_tables.py {input} > {output} 2>&1 | tee {log}"

# ======================================================================
# 25. Binning statistics (all bins)
# ======================================================================

rule binning_stats:
    input:
        metabat = "metabat/{sample}_bins",
        maxbin2 = "maxbin2/{sample}_bins",
        bin3c = "bin3c/{sample}/fasta",
        dastool = "dastool/{sample}/DASTool_bins" if config.get("run_dastool", True) else [],
        contigs = ASSEMBLY_CONTIGS_FILTERED
    output: "results/{sample}/binning_stats.txt"
    benchmark: "logs/binning_stats/{sample}.benchmark.txt"
    log: "logs/binning_stats/{sample}.log"
    conda: "envs/hicmag_py37.yaml"
    shell: "python scripts/binned_mags_stats.py {input.metabat} {input.maxbin2} {input.bin3c} {input.contigs} {output} 2>&1 | tee {log}"

# ======================================================================
# 26. MultiQC aggregation of all QC reports
# ======================================================================

rule multiqc_all:
    input:
        fastqc = expand("qc/fastqc/{sample}_WGS_R1_fastqc.zip", sample=SAMPLES),
        busco = expand("busco/{sample}/busco_results.tsv", sample=SAMPLES),
        quast = expand("results/{sample}/quast/report.html", sample=SAMPLES)
    output:
        html = "results/multiqc_report.html",
        done = touch("results/multiqc.done")
    benchmark: "logs/multiqc/benchmark.txt"
    log: "logs/multiqc/log.txt"
    conda: "envs/multiqc.yaml"
    shell:
        "multiqc -f "
        "-o results/multiqc "
        "-n multiqc_report "
        "qc/fastqc/ "
        "busco/ "
        "results/*/quast/ "
        "2>&1 | tee {log} && "
        "touch {output.done}"
