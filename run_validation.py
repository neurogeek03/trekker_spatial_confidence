#!/usr/bin/env python3
"""
Compute spatial coherence metrics and generate validation figures.

Tests whether confidence filtering improves the spatial clustering of
known cell types (DG granule, CA1/CA3 pyramidal, cortical L2/3, etc.).

Usage:
    python run_validation.py
    python run_validation.py --samples BC28,BC3 --fig-dir figures

Output:
    results/spatial_coherence_by_threshold.csv
    figures/validation_*.png (multiple figures)
"""
import argparse
import os
import numpy as np
import pandas as pd

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from spatial_confidence.io import load_merged_data
from spatial_confidence.validation import sweep_thresholds
from spatial_confidence.cell_types import VALIDATION_CELL_TYPES, HIPPO_CELL_TYPES
from spatial_confidence.plotting import (
    setup_plot_defaults, setup_spatial_axes, add_cell_count_label,
    plot_threshold_grid, build_legend_handles,
)
from spatial_confidence.config import CONFIDENCE_THRESHOLDS


def generate_coherence_metrics(sample_dfs, cell_types, results_dir):
    """Compute and save spatial coherence metrics across thresholds."""
    all_results = []
    for sid, df in sample_dfs.items():
        print(f"  Computing metrics for {sid}...")
        res = sweep_thresholds(df, cell_types, sample_id=sid)
        all_results.append(res)

    results_df = pd.concat(all_results, ignore_index=True)
    out_path = os.path.join(results_dir, 'spatial_coherence_by_threshold.csv')
    results_df.to_csv(out_path, index=False)
    print(f"  Saved: {out_path}")
    return results_df


def generate_hippocampal_spatial_maps(sample_dfs, fig_dir):
    """Generate hippocampal cell type spatial maps at different thresholds."""
    hippo_types = ['DG Granule', 'CA1 Pyramidal', 'CA3 Pyramidal']
    colors = {ct: HIPPO_CELL_TYPES[ct]['color'] for ct in hippo_types}
    thresholds = [0.0, 0.3, 0.5, 0.7]

    fig, axes = plt.subplots(len(sample_dfs), len(thresholds),
                             figsize=(28, 7 * len(sample_dfs)))
    if len(sample_dfs) == 1:
        axes = axes[np.newaxis, :]

    for row, (sid, df) in enumerate(sample_dfs.items()):
        pos = df[(df['n_clusters'] > 0) & df['subclass_name'].notna()]

        for col, thresh in enumerate(thresholds):
            ax = axes[row, col]
            filtered = pos[pos['confidence_score_penalized'] >= thresh]

            # Background
            other_mask = ~filtered['subclass_name'].fillna('').str.contains(
                '|'.join(info['pattern'] for info in HIPPO_CELL_TYPES.values()))
            ax.scatter(filtered[other_mask]['x_um'], filtered[other_mask]['y_um'],
                      s=0.3, alpha=0.08, color='#666666', rasterized=True)

            # Hippocampal types
            for ct_name in hippo_types:
                pattern = HIPPO_CELL_TYPES[ct_name]['pattern']
                cells = filtered[filtered['subclass_name'].fillna('').str.contains(pattern)]
                if len(cells) > 0:
                    ax.scatter(cells['x_um'], cells['y_um'], s=8, alpha=0.8,
                              color=colors[ct_name], rasterized=True,
                              label=f'{ct_name} ({len(cells):,})')

            setup_spatial_axes(ax)
            if row == 0:
                t = 'No filter' if thresh == 0.0 else f'Score ≥ {thresh}'
                ax.set_title(t, fontsize=20, fontweight='bold')
            if col == 0:
                ax.set_ylabel(f'{sid}\nY (µm)', fontsize=16)
            ax.legend(fontsize=9, loc='upper left', markerscale=2.5,
                     framealpha=0.8, facecolor='#1a1a2e', labelcolor='white')

    fig.suptitle('Hippocampal Cell Types at Different Confidence Thresholds',
                fontsize=24, fontweight='bold', y=0.98)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    path = os.path.join(fig_dir, 'validation_hippocampal_spatial_maps.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def generate_metrics_figure(results_df, sample_dfs, fig_dir):
    """Generate separation ratio and kNN distance vs threshold figure."""
    ct_colors = {ct: info['color'] for ct, info in VALIDATION_CELL_TYPES.items()}
    sample_ids = list(sample_dfs.keys())

    fig, axes = plt.subplots(2, len(sample_ids), figsize=(7 * len(sample_ids), 14))
    if len(sample_ids) == 1:
        axes = axes[:, np.newaxis]

    for col, sid in enumerate(sample_ids):
        # Separation ratio
        ax = axes[0, col]
        for ct_name, ct_info in VALIDATION_CELL_TYPES.items():
            sub = results_df[(results_df['cell_type'] == ct_name) & (results_df['sample'] == sid)]
            valid = sub.dropna(subset=['separation_ratio'])
            if len(valid) > 1:
                ax.plot(valid['threshold'], valid['separation_ratio'], 'o-',
                       color=ct_colors[ct_name], label=ct_name, linewidth=2.5, markersize=7)
        ax.axhline(1.0, color='white', ls=':', alpha=0.4)
        ax.set_xlabel('Confidence Threshold')
        ax.set_ylabel('Separation Ratio')
        ax.set_title(sid, fontsize=20, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.2)
        ax.set_facecolor('#f8f9fa')

        # kNN distance
        ax = axes[1, col]
        for ct_name, ct_info in VALIDATION_CELL_TYPES.items():
            sub = results_df[(results_df['cell_type'] == ct_name) & (results_df['sample'] == sid)]
            valid = sub.dropna(subset=['median_knn'])
            if len(valid) > 1:
                ax.plot(valid['threshold'], valid['median_knn'], 'o-',
                       color=ct_colors[ct_name], label=ct_name, linewidth=2.5, markersize=7)
        ax.set_xlabel('Confidence Threshold')
        ax.set_ylabel('Median Within-Type kNN (µm)')
        ax.set_title(sid, fontsize=20, fontweight='bold')
        ax.legend(fontsize=10)
        ax.grid(True, alpha=0.2)
        ax.set_facecolor('#f8f9fa')

    fig.suptitle('Spatial Coherence Improves with Confidence Filtering',
                fontsize=22, fontweight='bold', y=1.02)
    plt.tight_layout()
    path = os.path.join(fig_dir, 'validation_coherence_metrics.png')
    plt.savefig(path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {path}")


def main():
    parser = argparse.ArgumentParser(
        description='Compute spatial coherence metrics and generate validation figures.')
    parser.add_argument('--samples', default='BC28,BC3,BC13',
                        help='Comma-separated sample IDs')
    parser.add_argument('--results-dir', default='results')
    parser.add_argument('--fig-dir', default='figures')
    args = parser.parse_args()

    setup_plot_defaults()
    sample_ids = [s.strip() for s in args.samples.split(',')]
    os.makedirs(args.fig_dir, exist_ok=True)

    # Load data
    print("Loading merged data...")
    sample_dfs = {}
    for sid in sample_ids:
        sample_dfs[sid] = load_merged_data(args.results_dir, sid)
        n = len(sample_dfs[sid])
        n_ann = sample_dfs[sid]['subclass_name'].notna().sum()
        print(f"  {sid}: {n:,} cells, {n_ann:,} annotated")

    # Compute metrics
    print("\nComputing spatial coherence metrics...")
    results_df = generate_coherence_metrics(sample_dfs, VALIDATION_CELL_TYPES, args.results_dir)

    # Generate figures
    print("\nGenerating figures...")
    generate_hippocampal_spatial_maps(sample_dfs, args.fig_dir)
    generate_metrics_figure(results_df, sample_dfs, args.fig_dir)

    # Print summary
    print("\n" + "=" * 70)
    print("VALIDATION SUMMARY")
    print("=" * 70)
    for ct_name in ['DG Granule', 'CA1 Pyramidal', 'CA3 Pyramidal']:
        print(f"\n{ct_name}:")
        for sid in sample_ids:
            sub = results_df[(results_df['cell_type'] == ct_name) & (results_df['sample'] == sid)]
            t0 = sub[sub['threshold'] == 0.0]
            t5 = sub[sub['threshold'] == 0.5]
            if len(t0) > 0 and len(t5) > 0:
                r0 = t0['separation_ratio'].values[0]
                r5 = t5['separation_ratio'].values[0]
                if not np.isnan(r0) and not np.isnan(r5) and r0 > 0:
                    pct = (r5 - r0) / r0 * 100
                    print(f"  {sid}: sep {r0:.2f} → {r5:.2f} ({pct:+.0f}%)")

    print("\nDone!")


if __name__ == '__main__':
    main()
