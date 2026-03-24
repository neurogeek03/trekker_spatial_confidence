"""
Centralized configuration for the spatial confidence scoring pipeline.

All hardcoded constants live here. No other module should define these values.
To override for a specific run, pass parameters to the functions directly.
"""

# ── Trekker tile parameters ──────────────────────────────────────────────────

# Pixel-to-micron scale factor (Curio Trekker standard)
S = 0.647  # µm per pixel

# Maximum UMI count per bead-cell pair (Trekker's nUMI cap)
NUMI_MAX = 256

# ── DBSCAN parameters ────────────────────────────────────────────────────────

# Spatial radius for DBSCAN clustering (microns)
EPS = 50

# Default minPts for DBSCAN (Trekker uses 4)
MINPTS_DEFAULT = 4

# Cascade of minPts values to try: highest first, fall back to lower.
# minPts=4 is standard; 3 and 2 are used to rescue unpositioned cells.
MINPTS_CASCADE = [4, 3, 2]

# ── Confidence score parameters ──────────────────────────────────────────────

# Reference scale for cluster compactness normalization.
# Compactness = 1 / (1 + radius / COMPACTNESS_SCALE)
# A cluster with radius = COMPACTNESS_SCALE gets compactness = 0.5
COMPACTNESS_SCALE = 50  # µm

# Weights for each component in the composite confidence score.
# Components are rank-normalized to [0,1] before weighting.
COMPONENT_WEIGHTS = {
    'signal_fraction': 2.0,       # Most important: what fraction of UMI is signal
    'signal_n_beads': 1.5,        # More beads = better spatial triangulation
    'cluster_compactness': 1.5,   # Tighter cluster = more precise position
    'signal_gap_ratio_log': 1.0,  # Bigger gap between top clusters = less ambiguity
    'max_umi_enrichment': 1.0,    # Stronger peak signal = more confident
}

# Penalty factors for cells positioned at lower minPts thresholds.
# minPts=4 (standard) gets no penalty; lower values are less reliable.
MINPTS_PENALTIES = {4: 1.0, 3: 0.9, 2: 0.8}

# ── Validation parameters ────────────────────────────────────────────────────

# Number of nearest neighbors for spatial coherence metrics
K_NN = 10

# Standard confidence thresholds for sweep analyses
CONFIDENCE_THRESHOLDS = [0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7]

# ── Plotting defaults ────────────────────────────────────────────────────────

DARK_BG = '#0d1117'

PLOT_RCPARAMS = {
    'font.size': 14,
    'axes.titlesize': 20,
    'axes.labelsize': 16,
    'figure.facecolor': 'white',
    'figure.dpi': 100,
}

# ── Output columns ───────────────────────────────────────────────────────────

# Standard output columns for confidence score CSVs
OUTPUT_COLUMNS = [
    'cell_bc', 'confidence_score', 'confidence_score_penalized', 'minpts_used',
    'n_clusters', 'x_um', 'y_um',
    'signal_fraction', 'signal_n_beads', 'signal_umi',
    'cluster_radius_um', 'cluster_compactness',
    'noise_fraction', 'n_beads_total', 'total_umi', 'max_umi',
    'trekker_status',
]
