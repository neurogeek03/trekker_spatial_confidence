#!/usr/bin/env python3
"""
Generate specialized spatial visualization figures.

Usage:
    python run_visualization.py --type neurons    # all neuronal classes
    python run_visualization.py --type dcx        # DCX+/- DG granule cells
    python run_visualization.py --type all        # everything
"""
import argparse
import os
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from spatial_confidence.io import load_merged_data
from spatial_confidence.cell_types import (
    HIPPO_CELL_TYPES, build_neuronal_colormap,
)
from spatial_confidence.plotting import setup_plot_defaults, setup_spatial_axes, add_cell_count_label
from spatial_confidence.annotations import classify_neuronal


THRESHOLDS = [0.0, 0.3, 0.5, 0.7]


def generate_all_neurons(sample_dfs, fig_dir):
    """Generate all-neuronal-class spatial maps."""
    print("Generating all-neuron spatial maps...")

    # Collect all neuronal classes
    all_classes = set()
    for df in sample_dfs.values():
        neuronal = df[classify_neuronal(df)]
        all_classes.update(neuronal['class_name'].dropna().unique())
    all_classes = sorted(all_classes)
    class_colors = build_neuronal_colormap(all_classes)

    fig, axes = plt.subplots(len(sample_dfs), len(THRESHOLDS),
                             figsize=(32, 8 * len(sample_dfs)))
    if len(sample_dfs) == 1:
        axes = axes[np.newaxis, :]

    for row, (sid, df) in enumerate(sample_dfs.items()):
        pos = df[df['n_clusters'] > 0]
        ann = pos[pos['subclass_name'].notna()]
        neuronal_mask = classify_neuronal(ann)

        for col, thresh in enumerate(THRESHOLDS):
            ax = axes[row, col]
            filt = ann[ann['confidence_score_penalized'] >= thresh]
            n_all = len(pos[pos['confidence_score_penalized'] >= thresh])

            neurons = filt[neuronal_mask.reindex(filt.index, fill_value=False)]
            non_neurons = filt[~neuronal_mask.reindex(filt.index, fill_value=False)]

            ax.scatter(non_neurons['x_um'], non_neurons['y_um'],
                      s=0.3, alpha=0.06, color='#555555', rasterized=True)

            class_counts = neurons['class_name'].value_counts()
            for cls in class_counts.index[::-1]:
                cells = neurons[neurons['class_name'] == cls]
                ax.scatter(cells['x_um'], cells['y_um'], s=4, alpha=0.7,
                          color=class_colors.get(cls, '#fff'), rasterized=True)

            setup_spatial_axes(ax)
            add_cell_count_label(ax, f'{n_all:,} cells\n{len(neurons):,} neurons')

            if row == 0:
                t = 'No filter' if thresh == 0.0 else f'Score ≥ {thresh}'
                ax.set_title(t, fontsize=20, fontweight='bold')
            if col == 0:
                ax.set_ylabel(f'{sid}\nY (µm)', fontsize=16, fontweight='bold')

    fig.suptitle('All Neuronal Cell Types by Confidence Threshold',
                fontsize=26, fontweight='bold', y=1.01)
    plt.tight_layout(rect=[0, 0.02, 1, 0.99])
    path = os.path.join(fig_dir, 'all_neurons_spatial_combined.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def generate_dcx_dg(sample_dfs, fig_dir, dcx_dir='results'):
    """Generate DCX+/- DG granule cell figures."""
    print("Generating DCX+/- DG figures...")

    layers_config = [
        {'name': 'CA1 Pyramidal', 'query': lambda d: d['subclass_name'] == '016 CA1-ProS Glut',
         'color': '#2ecc71', 's': 6, 'alpha': 0.5, 'marker': 'o', 'zorder': 4},
        {'name': 'CA3 Pyramidal', 'query': lambda d: d['subclass_name'] == '017 CA3 Glut',
         'color': '#3498db', 's': 6, 'alpha': 0.5, 'marker': 'o', 'zorder': 4},
        {'name': 'DG DCX\u2212', 'query': lambda d: (d['subclass_name'] == '037 DG Glut') & (~d['DCX_positive']),
         'color': '#e74c3c', 's': 8, 'alpha': 0.6, 'marker': 'o', 'zorder': 5},
        {'name': 'DG DCX+', 'query': lambda d: (d['subclass_name'] == '037 DG Glut') & (d['DCX_positive']),
         'color': '#FFD700', 's': 120, 'alpha': 1.0, 'marker': '*', 'edgecolor': 'black', 'zorder': 12},
    ]

    # Merge DCX expression
    enriched_dfs = {}
    for sid, df in sample_dfs.items():
        dcx_path = os.path.join(dcx_dir, f'dcx_expression_{sid}.csv')
        if not os.path.exists(dcx_path):
            print(f"  WARNING: {dcx_path} not found. Run DCX extraction first.")
            continue
        dcx = pd.read_csv(dcx_path, index_col=0)
        merged = df.merge(dcx[['DCX_counts', 'DCX_positive']],
                         left_on='cell_bc', right_index=True, how='left')
        merged['DCX_positive'] = merged['DCX_positive'].fillna(False)
        enriched_dfs[sid] = merged

    if not enriched_dfs:
        print("  No DCX data available, skipping.")
        return

    # Full tile figure
    fig, axes = plt.subplots(len(enriched_dfs), len(THRESHOLDS),
                             figsize=(32, 8 * len(enriched_dfs)))
    if len(enriched_dfs) == 1:
        axes = axes[np.newaxis, :]

    for row, (sid, df) in enumerate(enriched_dfs.items()):
        pos = df[(df['n_clusters'] > 0) & df['subclass_name'].notna()]

        for col, thresh in enumerate(THRESHOLDS):
            ax = axes[row, col]
            filt = pos[pos['confidence_score_penalized'] >= thresh]
            n_all = len(df[(df['n_clusters'] > 0) & (df['confidence_score_penalized'] >= thresh)])

            # Background
            special_masks = [ly['query'](filt) for ly in layers_config]
            any_special = pd.concat(special_masks, axis=1).any(axis=1)
            other = filt[~any_special]
            ax.scatter(other['x_um'], other['y_um'], s=0.3, alpha=0.06,
                      color='#555555', rasterized=True)

            # Cell type layers
            counts = {}
            for ly in sorted(layers_config, key=lambda x: x['zorder']):
                cells = filt[ly['query'](filt)]
                counts[ly['name']] = len(cells)
                if len(cells) > 0:
                    kwargs = dict(s=ly['s'], alpha=ly['alpha'], color=ly['color'],
                                 marker=ly['marker'], zorder=ly['zorder'],
                                 rasterized=ly['s'] < 20)
                    if 'edgecolor' in ly:
                        kwargs['edgecolors'] = ly['edgecolor']
                        kwargs['linewidths'] = 1.0
                    ax.scatter(cells['x_um'], cells['y_um'], **kwargs)

            setup_spatial_axes(ax)
            add_cell_count_label(ax,
                f'{n_all:,} cells\nDCX+: {counts.get("DG DCX+",0)}  '
                f'DCX\u2212: {counts.get("DG DCX\u2212",0)}')

            if row == 0:
                t = 'No filter' if thresh == 0.0 else f'Score ≥ {thresh}'
                ax.set_title(t, fontsize=20, fontweight='bold')
            if col == 0:
                ax.set_ylabel(f'{sid}\nY (µm)', fontsize=16, fontweight='bold')

    # Legend
    legend_handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#e74c3c', markersize=9,
               label='DG Granule (DCX\u2212)', linestyle='None'),
        Line2D([0], [0], marker='*', color='w', markerfacecolor='#FFD700', markersize=18,
               markeredgecolor='black', markeredgewidth=1, label='DG Granule (DCX+)', linestyle='None'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#2ecc71', markersize=9,
               label='CA1 Pyramidal', linestyle='None'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor='#3498db', markersize=9,
               label='CA3 Pyramidal', linestyle='None'),
    ]
    fig.legend(legend_handles, [h.get_label() for h in legend_handles],
              loc='lower center', ncol=4, fontsize=14, bbox_to_anchor=(0.5, -0.03),
              frameon=True, facecolor='white', edgecolor='gray')

    fig.suptitle('DG Granule Cells by DCX Expression\nAcross Confidence Thresholds',
                fontsize=26, fontweight='bold', y=1.01)
    plt.tight_layout(rect=[0, 0.02, 1, 0.99])
    path = os.path.join(fig_dir, 'dcx_dg_spatial_clean.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def main():
    parser = argparse.ArgumentParser(
        description='Generate specialized spatial visualization figures.')
    parser.add_argument('--type', default='all', choices=['all', 'neurons', 'dcx'],
                        help='Figure type to generate')
    parser.add_argument('--samples', default='BC28,BC3,BC13')
    parser.add_argument('--results-dir', default='results')
    parser.add_argument('--fig-dir', default='figures')
    args = parser.parse_args()

    setup_plot_defaults()
    sample_ids = [s.strip() for s in args.samples.split(',')]
    os.makedirs(args.fig_dir, exist_ok=True)

    sample_dfs = {}
    for sid in sample_ids:
        sample_dfs[sid] = load_merged_data(args.results_dir, sid)
        print(f"  {sid}: {len(sample_dfs[sid]):,} cells")

    if args.type in ('all', 'neurons'):
        generate_all_neurons(sample_dfs, args.fig_dir)

    if args.type in ('all', 'dcx'):
        generate_dcx_dg(sample_dfs, args.fig_dir, dcx_dir=args.results_dir)

    print("\nDone!")


if __name__ == '__main__':
    main()
