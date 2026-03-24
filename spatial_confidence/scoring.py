"""
Core spatial confidence scoring functions.

Implements per-cell DBSCAN clustering, confidence component extraction,
rank-normalized composite scoring, and minPts penalty application.
"""
import numpy as np
import pandas as pd
from sklearn.cluster import DBSCAN
from scipy.stats import rankdata
import time

from .config import (
    EPS, MINPTS_DEFAULT, MINPTS_CASCADE, COMPACTNESS_SCALE,
    COMPONENT_WEIGHTS, MINPTS_PENALTIES, OUTPUT_COLUMNS,
)


def score_cell(grp, eps=EPS, minpts=MINPTS_DEFAULT,
               compactness_scale=COMPACTNESS_SCALE):
    """Run DBSCAN on a single cell's spatial barcodes and compute confidence metrics.

    Parameters
    ----------
    grp : pd.DataFrame
        Bead-level data for one cell. Must have columns: x_um, y_um, nUMI_collapsed.
    eps : float
        DBSCAN spatial radius in microns.
    minpts : int
        DBSCAN minimum cluster size.
    compactness_scale : float
        Reference radius for compactness normalization.

    Returns
    -------
    dict or None
        Dictionary of confidence metrics if DBSCAN finds ≥1 cluster, else None.
        Keys: n_clusters, x_um, y_um, signal_fraction, signal_n_beads, signal_umi,
        cluster_radius_um, cluster_compactness, signal_gap_ratio, max_umi_enrichment,
        noise_fraction, n_beads_total, total_umi, max_umi.
    """
    n_beads = len(grp)
    total_umi = grp['nUMI_collapsed'].sum()
    max_umi = grp['nUMI_collapsed'].max()

    if n_beads < minpts:
        return None

    X = grp[['x_um', 'y_um']].values
    umi = grp['nUMI_collapsed'].values

    db = DBSCAN(eps=eps, min_samples=minpts).fit(X)
    labels = db.labels_
    unique_labels = set(labels)
    n_clusters = len(unique_labels) - (1 if -1 in unique_labels else 0)

    if n_clusters == 0:
        return None

    # Find top cluster by UMI
    cluster_umi = {}
    cluster_n = {}
    for cl in unique_labels:
        if cl == -1:
            continue
        mask = labels == cl
        cluster_umi[cl] = umi[mask].sum()
        cluster_n[cl] = mask.sum()

    top_cluster = max(cluster_umi, key=cluster_umi.get)
    top_mask = labels == top_cluster
    top_umi = cluster_umi[top_cluster]
    top_n = cluster_n[top_cluster]
    top_X = X[top_mask]
    top_weights = umi[top_mask]

    # Signal fraction
    signal_fraction = top_umi / total_umi

    # UMI-weighted centroid
    centroid_x = np.average(top_X[:, 0], weights=top_weights)
    centroid_y = np.average(top_X[:, 1], weights=top_weights)

    # Cluster radius (weighted RMS distance from centroid)
    dists = np.sqrt((top_X[:, 0] - centroid_x)**2 + (top_X[:, 1] - centroid_y)**2)
    cluster_radius = np.sqrt(np.average(dists**2, weights=top_weights))

    # Compactness: 1/(1 + radius/scale), in [0, 1]
    compactness = 1.0 / (1.0 + cluster_radius / compactness_scale)

    # Signal gap: ratio of top to 2nd cluster (higher = less ambiguity)
    if n_clusters > 1:
        sorted_umis = sorted(cluster_umi.values(), reverse=True)
        signal_gap = sorted_umis[0] / max(sorted_umis[1], 1)
    else:
        signal_gap = top_umi  # single cluster = effectively infinite gap

    # Max UMI enrichment
    max_umi_enrichment = max_umi / max(total_umi / n_beads, 1e-10)

    return {
        'n_clusters': n_clusters,
        'x_um': centroid_x,
        'y_um': centroid_y,
        'signal_fraction': signal_fraction,
        'signal_n_beads': top_n,
        'signal_umi': top_umi,
        'cluster_radius_um': cluster_radius,
        'cluster_compactness': compactness,
        'signal_gap_ratio': signal_gap,
        'max_umi_enrichment': max_umi_enrichment,
        'noise_fraction': 1 - signal_fraction,
        'n_beads_total': n_beads,
        'total_umi': total_umi,
        'max_umi': max_umi,
    }


def score_all_cells(cell_groups, eps=EPS, minpts_cascade=None,
                    compactness_scale=COMPACTNESS_SCALE, verbose=True):
    """Score all cells using a minPts cascade strategy.

    For each cell, tries the highest minPts first and falls back to lower values.
    Cells that can't be clustered at any minPts get score=0.

    Parameters
    ----------
    cell_groups : dict
        Mapping of cell_bc → DataFrame (from load_trekker_data).
    eps : float
        DBSCAN spatial radius in microns.
    minpts_cascade : list of int, optional
        minPts values to try in order (default: [4, 3, 2]).
    compactness_scale : float
        Reference radius for compactness normalization.
    verbose : bool
        Print progress.

    Returns
    -------
    pd.DataFrame
        One row per cell with all confidence metric columns plus 'minpts_used'.
    """
    if minpts_cascade is None:
        minpts_cascade = list(MINPTS_CASCADE)

    results = []
    t0 = time.time()
    n_done = 0
    n_total = len(cell_groups)

    for cb in sorted(cell_groups.keys()):
        grp = cell_groups[cb]
        best = None
        best_minpts = 0

        for mp in minpts_cascade:
            res = score_cell(grp, eps, mp, compactness_scale)
            if res is not None:
                best = res
                best_minpts = mp
                break  # use the highest minPts that works

        if best is None:
            results.append({
                'cell_bc': cb,
                'minpts_used': 0,
                'n_clusters': 0,
                'x_um': np.nan, 'y_um': np.nan,
                'signal_fraction': 0, 'signal_n_beads': 0, 'signal_umi': 0,
                'cluster_radius_um': np.nan, 'cluster_compactness': 0,
                'signal_gap_ratio': 0, 'max_umi_enrichment': 0,
                'noise_fraction': 1, 'n_beads_total': len(grp),
                'total_umi': grp['nUMI_collapsed'].sum(),
                'max_umi': grp['nUMI_collapsed'].max(),
            })
        else:
            best['cell_bc'] = cb
            best['minpts_used'] = best_minpts
            results.append(best)

        n_done += 1
        if verbose and n_done % 5000 == 0:
            print(f"  Scored {n_done:,}/{n_total:,} cells...", flush=True)

    df_scores = pd.DataFrame(results)
    if verbose:
        elapsed = time.time() - t0
        n_positioned = (df_scores['n_clusters'] > 0).sum()
        print(f"  Done: {n_positioned:,}/{n_total:,} positioned in {elapsed:.0f}s", flush=True)

    return df_scores


def compute_composite_score(df_scores, component_weights=None,
                            minpts_penalties=None):
    """Compute the composite confidence score from individual components.

    Rank-normalizes each component within positioned cells, computes a weighted
    average, and applies a minPts penalty for cells found at lower thresholds.

    Parameters
    ----------
    df_scores : pd.DataFrame
        Output of score_all_cells. Modified in place.
    component_weights : dict, optional
        Weights for each component (default: COMPONENT_WEIGHTS).
    minpts_penalties : dict, optional
        Penalty multipliers by minPts value (default: MINPTS_PENALTIES).

    Returns
    -------
    pd.DataFrame
        Same DataFrame with added columns: confidence_score,
        confidence_score_penalized, minpts_penalty.
    """
    if component_weights is None:
        component_weights = dict(COMPONENT_WEIGHTS)
    if minpts_penalties is None:
        minpts_penalties = dict(MINPTS_PENALTIES)

    positioned = df_scores['n_clusters'] > 0
    df_scores['confidence_score'] = 0.0

    if positioned.sum() > 0:
        pos_idx = df_scores[positioned].index

        # Log-transform signal gap ratio (can be very large for single-cluster cells)
        df_scores['signal_gap_ratio_log'] = np.log1p(df_scores['signal_gap_ratio'])

        # Rank-normalize each component within positioned cells
        for comp in component_weights:
            vals = df_scores.loc[pos_idx, comp].values
            ranks = rankdata(vals, method='average')
            normed = (ranks - 1) / max(len(ranks) - 1, 1)
            df_scores.loc[pos_idx, f'{comp}_rank'] = normed

        # Weighted average
        total_weight = sum(component_weights.values())
        composite = np.zeros(len(pos_idx))
        for comp, weight in component_weights.items():
            composite += weight * df_scores.loc[pos_idx, f'{comp}_rank'].values
        composite /= total_weight
        df_scores.loc[pos_idx, 'confidence_score'] = composite

    # Apply minPts penalty
    df_scores['minpts_penalty'] = df_scores['minpts_used'].map(
        lambda mp: minpts_penalties.get(mp, 0.0)
    )
    df_scores['confidence_score_penalized'] = (
        df_scores['confidence_score'] * df_scores['minpts_penalty']
    )

    return df_scores


def add_trekker_status(df_scores, coords):
    """Add the original Trekker positioning status to scored cells.

    Parameters
    ----------
    df_scores : pd.DataFrame
        Output of compute_composite_score.
    coords : pd.DataFrame
        Trekker coords DataFrame with cell_bc and number_clusters columns.

    Returns
    -------
    pd.DataFrame
        Same DataFrame with added 'trekker_status' column.
    """
    trekker_clusters = coords.set_index('cell_bc')['number_clusters']
    df_scores['trekker_status'] = df_scores['cell_bc'].map(
        lambda cb: 'confident' if trekker_clusters.get(cb, 0) == 1 else
                   ('ambiguous' if trekker_clusters.get(cb, 0) > 1 else 'unpositioned')
    )
    return df_scores
