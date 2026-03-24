"""
Reusable plotting helpers for spatial confidence visualization.
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize

from .config import DARK_BG, PLOT_RCPARAMS, CONFIDENCE_THRESHOLDS


def setup_plot_defaults():
    """Apply standard matplotlib rcParams."""
    plt.rcParams.update(PLOT_RCPARAMS)


def setup_spatial_axes(ax, dark_bg=True, invert_y=True):
    """Configure axes for spatial scatter maps.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    dark_bg : bool
        Use dark background.
    invert_y : bool
        Invert Y axis (standard for tissue images).
    """
    if dark_bg:
        ax.set_facecolor(DARK_BG)
    ax.set_aspect('equal')
    if invert_y:
        ax.invert_yaxis()


def add_cell_count_label(ax, text, fontsize=11):
    """Add a white-on-black cell count annotation to axes.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    text : str
        Text to display (can include newlines).
    fontsize : int
    """
    ax.text(0.02, 0.02, text,
            transform=ax.transAxes, fontsize=fontsize, color='white',
            va='bottom', ha='left',
            bbox=dict(boxstyle='round,pad=0.3', facecolor='black', alpha=0.8))


def plot_spatial_scatter(ax, df, cell_type_layers, score_col='confidence_score_penalized',
                         threshold=0.0, show_background=True):
    """Plot a spatial map with multiple cell type layers.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
    df : pd.DataFrame
        Must have x_um, y_um, subclass_name, and score_col columns.
    cell_type_layers : list of dict
        Each dict has keys: name, query (callable), color, s, alpha, marker,
        edgecolor (optional), zorder.
    score_col : str
        Column to filter by threshold.
    threshold : float
        Minimum score to include.
    show_background : bool
        Plot non-typed cells as gray dots.

    Returns
    -------
    dict
        Cell counts per layer name.
    """
    filtered = df[df[score_col] >= threshold] if threshold > 0 else df

    # Identify all typed cells
    typed_masks = []
    for layer in cell_type_layers:
        if layer.get('query') is not None:
            typed_masks.append(layer['query'](filtered))

    if show_background and typed_masks:
        any_typed = pd.concat(typed_masks, axis=1).any(axis=1)
        other = filtered[~any_typed]
        ax.scatter(other['x_um'], other['y_um'], s=0.3, alpha=0.06,
                  color='#555555', rasterized=True)

    counts = {}
    for layer in cell_type_layers:
        if layer.get('query') is None:
            continue
        cells = filtered[layer['query'](filtered)]
        counts[layer['name']] = len(cells)
        if len(cells) > 0:
            kwargs = dict(
                s=layer.get('s', 6),
                alpha=layer.get('alpha', 0.6),
                color=layer['color'],
                marker=layer.get('marker', 'o'),
                zorder=layer.get('zorder', 5),
                rasterized=layer.get('s', 6) < 20,
            )
            ec = layer.get('edgecolor', 'none')
            if ec != 'none':
                kwargs['edgecolors'] = ec
                kwargs['linewidths'] = layer.get('linewidth', 1.0)
            ax.scatter(cells['x_um'], cells['y_um'], **kwargs)

    setup_spatial_axes(ax)
    return counts


def plot_threshold_grid(sample_dfs, cell_type_layers, thresholds=None,
                        fig_path=None, title=None, score_col='confidence_score_penalized',
                        figsize_per_panel=(8, 8)):
    """Generate the standard samples × thresholds grid figure.

    Parameters
    ----------
    sample_dfs : dict
        Mapping of {sample_label: DataFrame}.
    cell_type_layers : list of dict
        Cell type layer definitions (see plot_spatial_scatter).
    thresholds : list of float, optional
        Thresholds to show (default: [0.0, 0.3, 0.5, 0.7]).
    fig_path : str, optional
        Save path.
    title : str, optional
        Figure title.
    score_col : str
    figsize_per_panel : tuple

    Returns
    -------
    matplotlib.figure.Figure
    """
    if thresholds is None:
        thresholds = [0.0, 0.3, 0.5, 0.7]

    n_samples = len(sample_dfs)
    n_thresh = len(thresholds)
    figw = figsize_per_panel[0] * n_thresh
    figh = figsize_per_panel[1] * n_samples

    fig, axes = plt.subplots(n_samples, n_thresh, figsize=(figw, figh))
    if n_samples == 1:
        axes = axes[np.newaxis, :]

    for row, (label, df) in enumerate(sample_dfs.items()):
        pos = df[df['n_clusters'] > 0]
        ann = pos[pos['subclass_name'].notna()] if 'subclass_name' in pos.columns else pos

        for col, thresh in enumerate(thresholds):
            ax = axes[row, col]
            counts = plot_spatial_scatter(ax, ann, cell_type_layers,
                                         score_col=score_col, threshold=thresh)

            n_total = len(pos[pos[score_col] >= thresh])

            if row == 0:
                thresh_label = 'No filter' if thresh == 0.0 else f'Score ≥ {thresh}'
                ax.set_title(thresh_label, fontsize=20, fontweight='bold', pad=10)

            count_text = f'{n_total:,} cells'
            for name, n in counts.items():
                if n > 0:
                    count_text += f'\n{name}: {n}'
            add_cell_count_label(ax, count_text)

            if col == 0:
                ax.set_ylabel(f'{label}\nY (µm)', fontsize=16, fontweight='bold')
            if row == n_samples - 1:
                ax.set_xlabel('X (µm)', fontsize=14)

    if title:
        fig.suptitle(title, fontsize=26, fontweight='bold', y=1.01)

    plt.tight_layout(rect=[0, 0.02, 1, 0.99])

    if fig_path:
        plt.savefig(fig_path, dpi=150, bbox_inches='tight')
        plt.close()
        print(f"Saved: {fig_path}")

    return fig


def build_legend_handles(cell_type_layers):
    """Build matplotlib legend handles from cell type layer definitions.

    Parameters
    ----------
    cell_type_layers : list of dict

    Returns
    -------
    list of Line2D
    """
    handles = []
    for layer in cell_type_layers:
        if layer.get('query') is None:
            continue
        ms = 12 if layer.get('s', 6) >= 60 else 8
        ec = layer.get('edgecolor', 'none')
        h = Line2D([0], [0], marker=layer.get('marker', 'o'), color='w',
                   markerfacecolor=layer['color'], markersize=ms,
                   markeredgecolor=ec if ec != 'none' else layer['color'],
                   markeredgewidth=1 if ec != 'none' else 0,
                   label=layer.get('label', layer['name']), linestyle='None')
        handles.append(h)
    return handles
