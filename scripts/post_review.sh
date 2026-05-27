#!/bin/bash
# =============================================================================
# Wrapper: Run post-review pipeline for all samples or specific sample/MAG
# =============================================================================

set -e

show_usage() {
    cat <<EOF
Usage: $0 [OPTIONS]

Run post-Juicebox-review pipeline to regenerate corrected assemblies,
Hi‑C maps, statistics, and analysis.

Options:
  -s, --sample SAMPLE     Process specific sample only
  -m, --mag MAG           Process specific MAG only (requires -s)
  -c, --cores N           Number of cores (default: 32)
  -n, --dry-run           Show what would be done
  -h, --help              Show this help

Examples:
  $0                              # Run all samples
  $0 -s B32                      # Run sample B32 only
  $0 -s B32 -m bin.1.fa          # Run specific MAG
  $0 -n                          # Dry-run to see workflow
EOF
}

SAMPLE=""
MAG=""
CORES=32
DRY_RUN=false

while [[ $# -gt 0 ]]; do
    case $1 in
        -s|--sample) SAMPLE="$2"; shift 2 ;;
        -m|--mag) MAG="$2"; shift 2 ;;
        -c|--cores) CORES="$2"; shift 2 ;;
        -n|--dry-run) DRY_RUN=true; shift ;;
        -h|--help) show_usage; exit 0 ;;
        *) echo "Unknown option: $1"; show_usage; exit 1 ;;
    esac
done

# Build Snakemake command
SNAKEMAKE_CMD="snakemake --cores $CORES --use-conda -s Snakefile_post_review"

if [ "$DRY_RUN" = true ]; then
    SNAKEMAKE_CMD="$SNAKEMAKE_CMD -n"
fi

if [ -n "$SAMPLE" ] && [ -n "$MAG" ]; then
    TARGET="post_review/${SAMPLE}/${MAG}/post_review_complete.done"
elif [ -n "$SAMPLE" ]; then
    TARGET="post_review/${SAMPLE}/all_post_review.done"
else
    TARGET=""
fi

echo "========================================="
echo "Post-Review Pipeline"
echo "========================================="
echo "Sample: ${SAMPLE:-all}"
echo "MAG:    ${MAG:-all}"
echo "Cores:  $CORES"
echo "Dry-run: $DRY_RUN"
echo ""

if [ "$DRY_RUN" = true ]; then
    $SNAKEMAKE_CMD $TARGET
else
    $SNAKEMAKE_CMD $TARGET 2>&1 | tee post_review.log
    echo ""
    echo "Pipeline complete. Results in post_review/"
    echo "  - post_review/{sample}/{mag}/genome_corrected.fasta"
    echo "  - post_review/{sample}/{mag}/hic/inter_30_corrected.hic"
    echo "  - post_review/{sample}/{mag}/quast/report.html"
    echo "  - post_review/{sample}/{mag}/comparison_report.html"
fi
