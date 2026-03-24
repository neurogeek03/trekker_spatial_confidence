#!/usr/bin/env python3
"""
Merge cell type annotations with spatial confidence scores.

Loads annotations from an h5ad file and merges them with precomputed
confidence scores for each sample.

Usage:
    python run_annotations.py ../data/PCT_test_QC_merged_filtered_114914_mincells_10_in_2_samples_slide_tags.h5ad
    python run_annotations.py path/to/annotations.h5ad --samples BC28,BC3,BC13

Output:
    results/confidence_with_annotations_{SAMPLE_ID}.csv (one per sample)
"""
import argparse
import os

import matplotlib
matplotlib.use('Agg')

from spatial_confidence.io import load_annotations, load_confidence_scores
from spatial_confidence.annotations import merge_annotations


def main():
    parser = argparse.ArgumentParser(
        description='Merge cell type annotations with confidence scores.')
    parser.add_argument('h5ad_path', help='Path to annotated h5ad file')
    parser.add_argument('--samples', default='BC28,BC3,BC13',
                        help='Comma-separated sample IDs (default: BC28,BC3,BC13)')
    parser.add_argument('--results-dir', default='results',
                        help='Directory with confidence score CSVs (default: results)')
    args = parser.parse_args()

    sample_ids = [s.strip() for s in args.samples.split(',')]

    # Load annotations once
    ann_df = load_annotations(args.h5ad_path)

    # Merge for each sample
    for sid in sample_ids:
        print(f"\n--- {sid} ---")
        scores = load_confidence_scores(args.results_dir, sid, rescued=True)
        merged = merge_annotations(scores, ann_df, sample_id=sid)

        out_path = os.path.join(args.results_dir, f'confidence_with_annotations_{sid}.csv')
        merged.to_csv(out_path, index=False)
        print(f"  Saved: {out_path} ({len(merged):,} cells)")

    print("\nDone!")


if __name__ == '__main__':
    main()
