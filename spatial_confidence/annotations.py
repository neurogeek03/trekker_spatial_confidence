"""
Cell type annotation merging utilities.
"""
import pandas as pd


def merge_annotations(scores_df, ann_df, sample_id=None,
                      annotation_cols=None):
    """Merge cell type annotations with confidence scores.

    Parameters
    ----------
    scores_df : pd.DataFrame
        Confidence scores (from load_confidence_scores or score_all_cells).
    ann_df : pd.DataFrame
        Annotations (from load_annotations). Must have 'cell_bc' column
        and optionally 'sample' column.
    sample_id : str, optional
        If provided, filter ann_df to this sample before merging.
    annotation_cols : list of str, optional
        Columns to merge from ann_df. Default: class_name, subclass_label,
        subclass_name, subclass_bootstrapping_probability.

    Returns
    -------
    pd.DataFrame
        Merged DataFrame with all confidence columns plus annotations.
    """
    if annotation_cols is None:
        annotation_cols = ['class_name', 'subclass_label', 'subclass_name',
                           'subclass_bootstrapping_probability']

    ann = ann_df.copy()
    if sample_id is not None and 'sample' in ann.columns:
        ann = ann[ann['sample'] == sample_id]

    merge_cols = ['cell_bc'] + [c for c in annotation_cols if c in ann.columns]
    merged = scores_df.merge(ann[merge_cols], on='cell_bc', how='left')

    n_annotated = merged['subclass_name'].notna().sum() if 'subclass_name' in merged.columns else 0
    print(f"  Merged: {len(merged):,} cells, {n_annotated:,} annotated "
          f"({n_annotated/len(merged)*100:.1f}%)")
    return merged


def classify_neuronal(df, class_col='class_name', keywords=('Glut', 'GABA')):
    """Return a boolean Series indicating neuronal cells.

    Parameters
    ----------
    df : pd.DataFrame
        Must have the class_col column.
    class_col : str
        Column containing cell class names.
    keywords : tuple of str
        Keywords to match for neuronal classes.

    Returns
    -------
    pd.Series of bool
    """
    return df[class_col].fillna('').apply(
        lambda x: any(k in str(x) for k in keywords)
    )
