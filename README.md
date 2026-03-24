# trekker_spatial_confidence

Continuous spatial positioning confidence scoring for [Curio Trekker](https://curiobioscience.com/) / Slide-Tags spatial transcriptomics data.

## What this does

The Curio Trekker pipeline classifies cells as "confidently positioned" (1 DBSCAN cluster), "ambiguous" (2+ clusters), or "unpositioned" (0 clusters). This tool computes a **continuous confidence score (0-1)** for each cell, enabling quality-based filtering that dramatically improves spatial analysis.

**Key result**: Filtering to score >= 0.5 produces **3-8x improvement** in spatial clustering metrics for hippocampal cell types (DG granule, CA1, CA3 pyramidal cells).

See [REPORT.md](REPORT.md) for full methodology, validation, and recommendations.

## Quick Start

```bash
# Score a sample (~3 min on a laptop)
python run_scoring.py BC28 /path/to/trekker/output

# Merge cell type annotations
python run_annotations.py /path/to/annotations.h5ad --samples BC28

# Generate validation figures
python run_validation.py --samples BC28
```

## Requirements

```
numpy
pandas
scikit-learn
scipy
matplotlib
anndata
```

## Required Input Files (per sample)

From Trekker output:
- `df_whitelist_{SAMPLE_ID}.txt` (spatial barcode UMI counts)
- `matching_result_{SAMPLE_ID}.csv` (barcode-to-bead coordinate mapping)
- `coords_{SAMPLE_ID}.txt` (Trekker DBSCAN results)

## Confidence Score Components

| Component | Weight | Measures |
|-----------|--------|----------|
| Signal fraction | 2.0 | % of UMI in top DBSCAN cluster |
| Signal bead count | 1.5 | Number of beads in top cluster |
| Cluster compactness | 1.5 | Inverse spatial radius of cluster |
| Signal gap ratio | 1.0 | Dominance of top vs 2nd cluster |
| Max UMI enrichment | 1.0 | Peak barcode signal strength |

## Recommended Thresholds

| Threshold | Retention | Use case |
|-----------|----------|----------|
| >= 0.3 | ~73% | Exploratory analyses |
| >= 0.5 | ~42% | **Spatial precision analyses** |
| >= 0.7 | ~22% | Maximum spatial quality |

## Project Structure

```
spatial_confidence/     # Python package
    config.py           # All constants and parameters
    io.py               # Data loading
    scoring.py          # DBSCAN + confidence scoring
    annotations.py      # Cell type annotation merging
    validation.py       # Spatial coherence metrics
    plotting.py         # Visualization helpers
    cell_types.py       # Cell type definitions
run_scoring.py          # Score a sample
run_annotations.py      # Merge annotations
run_validation.py       # Validation metrics + figures
run_visualization.py    # Specialized figures
```
