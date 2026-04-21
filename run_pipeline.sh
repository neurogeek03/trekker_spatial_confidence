#!/usr/bin/env bash
set -euo pipefail

# Run the full trekker_spatial_confidence pipeline across all samples.
# Must be run from trekker_spatial_confidence/ (the directory containing this script).
#
# Usage:
#   ./run_pipeline.sh                      # all 4 stages
#   ./run_pipeline.sh --stage scoring
#   ./run_pipeline.sh --stage annotations
#   ./run_pipeline.sh --stage validation
#   ./run_pipeline.sh --stage visualization

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/config/pipeline.conf"

# Parse --stage argument (default: all)
STAGE="all"
while [[ $# -gt 0 ]]; do
    case "$1" in
        --stage) STAGE="$2"; shift 2 ;;
        *) echo "Unknown argument: $1" >&2; exit 1 ;;
    esac
done

# Run from the script's own directory
cd "$SCRIPT_DIR"

# Activate venv
source "$VENV/bin/activate"

# Create output directories
mkdir -p "$RESULTS_DIR" "$FIGURES_DIR" "$LOG_DIR"

# Log file
LOGFILE="$LOG_DIR/pipeline_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee "$LOGFILE") 2>&1

echo "============================================================"
echo "Trekker Spatial Confidence Pipeline"
echo "Started: $(date)"
echo "Stage:   $STAGE"
echo "Samples: ${SAMPLES[*]}"
echo "============================================================"

# Comma-separated sample list (for scripts that take --samples)
SAMPLES_CSV=$(IFS=,; echo "${SAMPLES[*]}")

# ── Stage 1: Scoring ─────────────────────────────────────────────────────────
run_scoring() {
    echo ""
    echo "--- STAGE 1: Scoring ---"
    for SAMPLE in "${SAMPLES[@]}"; do
        echo "[$(date +%H:%M:%S)] Scoring $SAMPLE..."
        python run_scoring.py "$SAMPLE" "$DATA_DIR/$SAMPLE" \
            --output-dir "$RESULTS_DIR" \
            --fig-dir "$FIGURES_DIR"
    done
    echo "[$(date +%H:%M:%S)] Scoring complete."
}

# ── Stage 2: Annotations ─────────────────────────────────────────────────────
run_annotations() {
    echo ""
    echo "--- STAGE 2: Annotations ---"
    echo "[$(date +%H:%M:%S)] Merging annotations for: $SAMPLES_CSV"
    python run_annotations.py "$ANNOTATIONS_H5AD" \
        --samples "$SAMPLES_CSV" \
        --results-dir "$RESULTS_DIR"
    echo "[$(date +%H:%M:%S)] Annotations complete."
}

# ── Stage 3: Validation ──────────────────────────────────────────────────────
run_validation() {
    echo ""
    echo "--- STAGE 3: Validation ---"
    echo "[$(date +%H:%M:%S)] Running validation for: $SAMPLES_CSV"
    python run_validation.py \
        --samples "$SAMPLES_CSV" \
        --results-dir "$RESULTS_DIR" \
        --fig-dir "$FIGURES_DIR"
    echo "[$(date +%H:%M:%S)] Validation complete."
}

# ── Stage 4: Visualization ───────────────────────────────────────────────────
run_visualization() {
    echo ""
    echo "--- STAGE 4: Visualization ---"
    echo "[$(date +%H:%M:%S)] Running visualization for: $SAMPLES_CSV"
    python run_visualization.py \
        --samples "$SAMPLES_CSV" \
        --results-dir "$RESULTS_DIR" \
        --fig-dir "$FIGURES_DIR"
    echo "[$(date +%H:%M:%S)] Visualization complete."
}

# ── Dispatch ─────────────────────────────────────────────────────────────────
case "$STAGE" in
    all)
        run_scoring
        run_annotations
        run_validation
        run_visualization
        ;;
    scoring)       run_scoring ;;
    annotations)   run_annotations ;;
    validation)    run_validation ;;
    visualization) run_visualization ;;
    *)
        echo "Unknown stage: $STAGE" >&2
        echo "Valid stages: all, scoring, annotations, validation, visualization" >&2
        exit 1
        ;;
esac

echo ""
echo "============================================================"
echo "Done: $(date)"
echo "Log: $LOGFILE"
echo "============================================================"
