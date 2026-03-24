"""
Spatial coherence validation metrics.

Computes separation ratio and within-type kNN distance to quantify
how well cell types cluster spatially at different confidence thresholds.
"""
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors

from .config import K_NN, CONFIDENCE_THRESHOLDS


def compute_spatial_metrics(xy, k=K_NN):
    """Compute within-type spatial metrics for a set of cell positions.

    Parameters
    ----------
    xy : np.ndarray, shape (n, 2)
        X, Y coordinates.
    k : int
        Number of nearest neighbors.

    Returns
    -------
    dict
        Keys: median_knn, mean_knn, spatial_extent, n_cells.
        Returns NaN values if too few cells.
    """
    if len(xy) < k + 1:
        return {'median_knn': np.nan, 'mean_knn': np.nan,
                'spatial_extent': np.nan, 'n_cells': len(xy)}

    nn = NearestNeighbors(n_neighbors=k + 1)
    nn.fit(xy)
    dists, _ = nn.kneighbors(xy)
    knn_dists = dists[:, 1:]  # exclude self

    per_cell_mean_knn = knn_dists.mean(axis=1)

    # Spatial extent: IQR-based (robust to outliers)
    x_iqr = np.percentile(xy[:, 0], 75) - np.percentile(xy[:, 0], 25)
    y_iqr = np.percentile(xy[:, 1], 75) - np.percentile(xy[:, 1], 25)

    return {
        'median_knn': np.median(per_cell_mean_knn),
        'mean_knn': np.mean(per_cell_mean_knn),
        'spatial_extent': x_iqr * y_iqr,
        'n_cells': len(xy),
    }


def compute_separation_ratio(xy_type, xy_other, k=K_NN):
    """Compute the spatial separation of a cell type from other cells.

    Ratio = (mean distance to k nearest OTHER-type cells) /
            (mean distance to k nearest SAME-type cells).
    Values > 1 indicate spatial clustering; higher is better.

    Parameters
    ----------
    xy_type : np.ndarray, shape (n, 2)
        Positions of cells of the target type.
    xy_other : np.ndarray, shape (m, 2)
        Positions of all other cells.
    k : int
        Number of nearest neighbors.

    Returns
    -------
    float
        Median separation ratio across cells of the target type.
        Returns NaN if insufficient cells.
    """
    if len(xy_type) < k + 1 or len(xy_other) < k:
        return np.nan

    # Distance to same-type neighbors
    nn_same = NearestNeighbors(n_neighbors=min(k + 1, len(xy_type)))
    nn_same.fit(xy_type)
    dists_same, _ = nn_same.kneighbors(xy_type)
    mean_same = dists_same[:, 1:].mean(axis=1)

    # Distance to other-type neighbors
    nn_other = NearestNeighbors(n_neighbors=min(k, len(xy_other)))
    nn_other.fit(xy_other)
    dists_other, _ = nn_other.kneighbors(xy_type)
    mean_other = dists_other.mean(axis=1)

    ratio = mean_other / (mean_same + 1e-6)
    return np.median(ratio)


def sweep_thresholds(df, cell_types, sample_id='',
                     thresholds=None, k=K_NN,
                     score_col='confidence_score_penalized'):
    """Run spatial coherence analysis across confidence thresholds.

    For each threshold × cell type, computes separation ratio and kNN distance.

    Parameters
    ----------
    df : pd.DataFrame
        Merged confidence + annotation data. Must have columns: n_clusters,
        x_um, y_um, subclass_name, and score_col.
    cell_types : dict
        Mapping of display name → dict with 'pattern' key (matched against
        subclass_name via str.contains).
    sample_id : str
        Label for the sample column in output.
    thresholds : list of float, optional
        Confidence thresholds to sweep (default: CONFIDENCE_THRESHOLDS).
    k : int
        Number of nearest neighbors.
    score_col : str
        Column name for the confidence score.

    Returns
    -------
    pd.DataFrame
        Columns: sample, threshold, cell_type, n_cells, median_knn,
        mean_knn, spatial_extent, separation_ratio.
    """
    if thresholds is None:
        thresholds = list(CONFIDENCE_THRESHOLDS)

    positioned = df[(df['n_clusters'] > 0) & df['x_um'].notna() & df['subclass_name'].notna()]
    results = []

    for thresh in thresholds:
        filtered = positioned[positioned[score_col] >= thresh]

        for ct_name, ct_info in cell_types.items():
            ct_mask = filtered['subclass_name'].fillna('').str.contains(ct_info['pattern'])
            ct_cells = filtered[ct_mask]
            other_cells = filtered[~ct_mask]

            if len(ct_cells) < k + 1:
                results.append({
                    'sample': sample_id, 'threshold': thresh, 'cell_type': ct_name,
                    'n_cells': len(ct_cells), 'median_knn': np.nan,
                    'separation_ratio': np.nan, 'spatial_extent': np.nan,
                })
                continue

            xy_type = ct_cells[['x_um', 'y_um']].values
            xy_other = other_cells[['x_um', 'y_um']].values

            metrics = compute_spatial_metrics(xy_type, k=k)
            sep = compute_separation_ratio(xy_type, xy_other, k=k)

            results.append({
                'sample': sample_id, 'threshold': thresh, 'cell_type': ct_name,
                'n_cells': len(ct_cells),
                'median_knn': metrics['median_knn'],
                'mean_knn': metrics['mean_knn'],
                'spatial_extent': metrics['spatial_extent'],
                'separation_ratio': sep,
            })

    return pd.DataFrame(results)
