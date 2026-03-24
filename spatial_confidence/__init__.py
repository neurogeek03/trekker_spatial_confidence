"""
Spatial confidence scoring for Curio Trekker / Slide-Tags data.

Computes a continuous spatial positioning confidence score for each cell,
enabling quality-based filtering for downstream spatial analyses.

Usage:
    from spatial_confidence import score_sample, merge_annotations

See run_scoring.py for full pipeline entry point.
"""
from .scoring import score_cell, score_all_cells, compute_composite_score, add_trekker_status
from .io import load_trekker_data, load_confidence_scores, load_annotations
from .annotations import merge_annotations
