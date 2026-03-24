#!/usr/bin/env python3
"""
Score a Trekker sample with spatial confidence metrics.

This is the main entry point for computing per-cell spatial positioning
confidence scores. It replaces the previous 05_confidence_score.py and
09_score_unpositioned.py scripts.

Usage:
    python run_scoring.py BC28 ../data
    python run_scoring.py BC13 ../data/BC13/BC13_bad_sample
    python run_scoring.py BC28 ../data --minpts-cascade 4  # standard only, no rescue

Output:
    results/confidence_scores_rescued_{SAMPLE_ID}.csv
    figures/rescue_unpositioned_{SAMPLE_ID}.png
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from spatial_confidence.io import load_trekker_data
from spatial_confidence.scoring import score_all_cells, compute_composite_score, add_trekker_status
from spatial_confidence.config import OUTPUT_COLUMNS


def print_summary(df_scores):
    """Print a summary of scoring results."""
    print("\n" + "=" * 70)
    print("SCORING SUMMARY")
    print("=" * 70)

    for mp in sorted(df_scores['minpts_used'].unique(), reverse=True):
        n = (df_scores['minpts_used'] == mp).sum()
        label = f'minPts={mp}' if mp > 0 else 'No cluster'
        print(f"  {label:15s}: {n:,} cells")

    print(f"\nBy Trekker status:")
    for status in ['confident', 'ambiguous', 'unpositioned']:
        sub = df_scores[df_scores['trekker_status'] == status]
        scored = sub[sub['confidence_score'] > 0]
        if len(scored) > 0:
            print(f"  {status:15s}: {len(sub):,} total, {len(scored):,} scored, "
                  f"median={scored['confidence_score_penalized'].median():.3f}")
        else:
            print(f"  {status:15s}: {len(sub):,} total, 0 scored")

    print(f"\nThreshold summary:")
    for thresh in [0.3, 0.4, 0.5, 0.6, 0.7]:
        above = (df_scores['confidence_score_penalized'] >= thresh).sum()
        pct = above / len(df_scores) * 100
        print(f"  score >= {thresh}: {above:,} cells ({pct:.1f}%)")


def save_diagnostic_figure(df_scores, sample_id, fig_dir):
    """Save a 6-panel diagnostic figure."""
    fig, axes = plt.subplots(2, 3, figsize=(20, 12))

    colors_status = {'confident': 'green', 'ambiguous': 'orange', 'unpositioned': 'red'}
    colors_mp = {4: 'green', 3: 'orange', 2: 'red'}

    # 1. Score by Trekker status
    ax = axes[0, 0]
    for status in ['confident', 'ambiguous', 'unpositioned']:
        sub = df_scores[(df_scores['trekker_status'] == status) &
                        (df_scores['confidence_score_penalized'] > 0)]
        if len(sub) > 0:
            ax.hist(sub['confidence_score_penalized'], bins=60, alpha=0.4,
                    color=colors_status[status], label=f'{status} ({len(sub):,})', density=True)
    ax.set_xlabel('Penalized Confidence Score')
    ax.set_ylabel('Density')
    ax.set_title('Score by Trekker status')
    ax.legend(fontsize=9)

    # 2. Score by minPts
    ax = axes[0, 1]
    for mp in [4, 3, 2]:
        sub = df_scores[(df_scores['minpts_used'] == mp) &
                        (df_scores['confidence_score_penalized'] > 0)]
        if len(sub) > 0:
            ax.hist(sub['confidence_score_penalized'], bins=60, alpha=0.4,
                    color=colors_mp[mp], label=f'minPts={mp} ({len(sub):,})', density=True)
    ax.set_xlabel('Penalized Confidence Score')
    ax.set_ylabel('Density')
    ax.set_title('Score by minPts used')
    ax.legend(fontsize=9)

    # 3. Signal fraction by minPts
    import numpy as np
    ax = axes[0, 2]
    for mp in [4, 3, 2]:
        sub = df_scores[df_scores['minpts_used'] == mp]
        if len(sub) > 0:
            ax.hist(np.log10(sub['signal_fraction'].clip(lower=1e-6)), bins=60, alpha=0.4,
                    color=colors_mp[mp], label=f'minPts={mp}', density=True)
    ax.set_xlabel('log10(Signal Fraction)')
    ax.set_ylabel('Density')
    ax.set_title('Signal fraction by minPts')
    ax.legend(fontsize=9)

    # 4-6. Spatial maps
    for idx, (mask_label, mask_fn, title_fn) in enumerate([
        ('minPts=4', lambda d: (d['minpts_used'] == 4) & d['x_um'].notna(),
         lambda n: f'minPts=4 only ({n:,})'),
        ('All scored', lambda d: (d['n_clusters'] > 0) & d['x_um'].notna(),
         lambda n: f'All scored ({n:,})'),
        ('Rescued', lambda d: (d['trekker_status'] == 'unpositioned') &
                              (d['n_clusters'] > 0) & d['x_um'].notna(),
         lambda n: f'Rescued ({n:,})'),
    ]):
        ax = axes[1, idx]
        ax.set_facecolor('black')
        sub = df_scores[mask_fn(df_scores)].sort_values('confidence_score_penalized')
        if len(sub) > 0:
            ax.scatter(sub['x_um'], sub['y_um'], c=sub['confidence_score_penalized'],
                      s=2, alpha=0.5, cmap='inferno', vmin=0, vmax=1, rasterized=True)
        ax.set_title(title_fn(len(sub)), fontsize=14)
        ax.set_aspect('equal')

    plt.suptitle(f'Spatial Confidence Score: {sample_id}', fontsize=22, y=1.02)
    plt.tight_layout()
    path = os.path.join(fig_dir, f'rescue_unpositioned_{sample_id}.png')
    plt.savefig(path, dpi=200, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"Saved: {path}")


def main():
    parser = argparse.ArgumentParser(
        description='Compute spatial confidence scores for a Trekker sample.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  python run_scoring.py BC28 ../data
  python run_scoring.py BC13 ../data/BC13/BC13_bad_sample
  python run_scoring.py BC28 ../data --minpts-cascade 4    # no rescue
        """)
    parser.add_argument('sample_id', help='Sample ID (e.g. BC28)')
    parser.add_argument('data_dir', help='Path to Trekker output directory')
    parser.add_argument('--minpts-cascade', default='4,3,2',
                        help='Comma-separated minPts values to try (default: 4,3,2)')
    parser.add_argument('--output-dir', default='results',
                        help='Output directory for CSVs (default: results)')
    parser.add_argument('--fig-dir', default='figures',
                        help='Output directory for figures (default: figures)')
    args = parser.parse_args()

    cascade = [int(x) for x in args.minpts_cascade.split(',')]
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs(args.fig_dir, exist_ok=True)

    # 1. Load data
    coords, cell_groups = load_trekker_data(args.data_dir, args.sample_id)

    # 2. Score all cells
    df_scores = score_all_cells(cell_groups, minpts_cascade=cascade)

    # 3. Composite score + penalties
    df_scores = compute_composite_score(df_scores)

    # 4. Merge Trekker status
    df_scores = add_trekker_status(df_scores, coords)

    # 5. Save
    out_cols = [c for c in OUTPUT_COLUMNS if c in df_scores.columns]
    out_path = os.path.join(args.output_dir, f'confidence_scores_rescued_{args.sample_id}.csv')
    df_scores[out_cols].to_csv(out_path, index=False)
    print(f"Saved: {out_path}")

    # 6. Summary + figure
    print_summary(df_scores)
    save_diagnostic_figure(df_scores, args.sample_id, args.fig_dir)


if __name__ == '__main__':
    main()
