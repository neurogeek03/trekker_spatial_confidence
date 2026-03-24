"""
Data loading utilities for the spatial confidence scoring pipeline.

Handles loading and preprocessing of Trekker output files (coords, df_whitelist,
matching_result) and CellBender/annotation h5ad files.
"""
import os
import time
import numpy as np
import pandas as pd
from .config import S, NUMI_MAX


def find_file(data_dir, sample_id, filename_pattern, subdirs=('', 'misc')):
    """Find a file in data_dir or its subdirectories.

    Searches for files matching filename_pattern.format(sample_id=sample_id)
    in data_dir and then in data_dir/subdir for each subdir.

    Parameters
    ----------
    data_dir : str
        Base directory to search.
    sample_id : str
        Sample identifier to substitute into filename_pattern.
    filename_pattern : str
        Filename with {sample_id} placeholder, e.g. 'coords_{sample_id}.txt'
    subdirs : tuple of str
        Subdirectories to search ('' means data_dir itself).

    Returns
    -------
    str
        Full path to the found file.

    Raises
    ------
    FileNotFoundError
        If the file is not found in any searched location.
    """
    filename = filename_pattern.format(sample_id=sample_id)
    searched = []
    for sub in subdirs:
        path = os.path.join(data_dir, sub, filename) if sub else os.path.join(data_dir, filename)
        searched.append(path)
        if os.path.exists(path):
            return path
    raise FileNotFoundError(
        f"Cannot find {filename} in: {searched}"
    )


def load_trekker_data(data_dir, sample_id, s=S, numi_max=NUMI_MAX, verbose=True):
    """Load and preprocess Trekker output files for confidence scoring.

    Performs the Trekker-faithful preprocessing:
    1. Load coords (DBSCAN results from Trekker)
    2. Load df_whitelist (spatial barcode UMI per cell)
    3. Load matching_result (spatial barcode → bead coordinates)
    4. Merge, collapse by bead, convert pixel→micron, filter by nUMI cap

    Parameters
    ----------
    data_dir : str
        Directory containing Trekker output files.
    sample_id : str
        Sample identifier (e.g. 'BC28').
    s : float
        Pixel-to-micron scale factor (default: 0.647).
    numi_max : int
        Maximum UMI count per bead-cell pair (default: 256).
    verbose : bool
        Print progress messages.

    Returns
    -------
    coords : pd.DataFrame
        Trekker coordinates with columns: cell_bc, x_um, y_um, number_clusters, etc.
    cell_groups : dict
        Mapping of cell_bc → DataFrame with columns: matched_beadbarcode,
        nUMI_collapsed, x_um, y_um. Ready for DBSCAN.
    """
    if verbose:
        print(f"Loading Trekker data for {sample_id}...", flush=True)
    t0 = time.time()

    # 1. Coords
    coords_path = find_file(data_dir, sample_id, 'coords_{sample_id}.txt')
    coords = pd.read_csv(coords_path, sep=r'\s+')
    n_total = len(coords)
    if verbose:
        n_conf = (coords['number_clusters'] == 1).sum()
        n_unpos = (coords['number_clusters'] == 0).sum()
        print(f"  Coords: {n_total:,} cells — {n_conf:,} confident, {n_unpos:,} unpositioned")

    # 2. Whitelist
    wl_path = find_file(data_dir, sample_id, 'df_whitelist_{sample_id}.txt')
    df_wl = pd.read_csv(wl_path)
    if verbose:
        print(f"  Whitelist: {len(df_wl):,} rows")

    # 3. Matching result
    match_path = find_file(data_dir, sample_id, 'matching_result_{sample_id}.csv')
    matching = pd.read_csv(match_path)
    matching['x_um'] = matching['xcoord'] * s
    matching['y_um'] = matching['ycoord'] * s
    sb_matched = matching[['matched_beadbarcode', 'Illumina_barcode', 'x_um', 'y_um']].copy()

    # 4. Trekker-style merge + collapse
    cb_sb = sb_matched.merge(df_wl, left_on='Illumina_barcode', right_on='SB', how='inner').drop_duplicates()
    collapsed = cb_sb.groupby(['matched_beadbarcode', 'CB']).agg(
        nUMI_collapsed=('nUMI', 'sum')
    ).reset_index()
    bead_coords_df = cb_sb[['matched_beadbarcode', 'x_um', 'y_um']].drop_duplicates(
        subset='matched_beadbarcode', keep='first')
    collapsed = collapsed.merge(bead_coords_df, on='matched_beadbarcode', how='left')
    collapsed = collapsed[collapsed['nUMI_collapsed'] < numi_max]

    # 5. Group by cell
    cell_groups = {cb: grp for cb, grp in collapsed.groupby('CB')}

    if verbose:
        print(f"  Processed {len(cell_groups):,} cells in {time.time()-t0:.1f}s", flush=True)

    return coords, cell_groups


def load_confidence_scores(results_dir, sample_id, rescued=True):
    """Load precomputed confidence scores from CSV.

    Parameters
    ----------
    results_dir : str
        Directory containing result CSVs.
    sample_id : str
        Sample identifier.
    rescued : bool
        If True, load the rescued version (with minPts cascade); otherwise load
        the standard version.

    Returns
    -------
    pd.DataFrame
    """
    prefix = 'confidence_scores_rescued' if rescued else 'confidence_scores'
    path = os.path.join(results_dir, f'{prefix}_{sample_id}.csv')
    return pd.read_csv(path)


def load_annotations(h5ad_path, annotation_cols=None):
    """Load cell type annotations from an h5ad file.

    Strips the '-1' suffix from barcodes to match Trekker barcode format.

    Parameters
    ----------
    h5ad_path : str
        Path to the annotated h5ad file.
    annotation_cols : list of str, optional
        Columns to extract from obs. Default: sample, class_name, subclass_label,
        subclass_name, subclass_bootstrapping_probability.

    Returns
    -------
    pd.DataFrame
        DataFrame with cell_bc column and annotation columns.
    """
    import anndata as ad

    if annotation_cols is None:
        annotation_cols = ['sample', 'class_name', 'subclass_label', 'subclass_name',
                           'subclass_bootstrapping_probability']

    print(f"Loading annotations from {h5ad_path}...", flush=True)
    adata = ad.read_h5ad(h5ad_path, backed='r')
    ann_df = adata.obs[annotation_cols].copy()
    ann_df.index = ann_df.index.str.replace('-1$', '', regex=True)
    ann_df = ann_df.reset_index(names='cell_bc')
    print(f"  {len(ann_df):,} cells loaded")
    return ann_df


def load_merged_data(results_dir, sample_id):
    """Load confidence scores merged with annotations.

    Parameters
    ----------
    results_dir : str
        Directory containing result CSVs.
    sample_id : str
        Sample identifier.

    Returns
    -------
    pd.DataFrame
    """
    path = os.path.join(results_dir, f'confidence_with_annotations_{sample_id}.csv')
    return pd.read_csv(path)
