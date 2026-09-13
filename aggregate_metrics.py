#!/usr/bin/env python3
# scripts/aggregate_metrics.py
# Aggregate per-sample metrics and generate summary tables Safely

import os
import sys
import pandas as pd
import numpy as np


def main():
    if len(sys.argv) != 2:
        print("Usage: aggregate_metrics.py <results/validation/all_metrics.tsv>")
        sys.exit(1)

    metrics_file = sys.argv[1]
    ROOT_DIR = os.path.dirname(metrics_file)

    # --- Load metrics table ------------------------------------------------
    df = pd.read_csv(metrics_file, sep='\t', na_values=['NA'])
    df = df.dropna(subset=['Tool'])

    # --- Coerce all numeric columns ----------------------------------------
    numeric_cols = [
        'Ref_Plasmid_Count', 'Pred_Plasmid_Count',
        'TP', 'FP', 'FN',
        'Circular_TP', 'Circularization_Pct',
        'Accuracy', 'Sensitivity',
        'Contig_Precision', 'BP_Precision', 'F1_Score',
        'Mean_Identity', 'Mean_Completeness', 'Contamination_Pct'
    ]
    for col in numeric_cols:
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    # ===========================================================================
    # SAFE-GUARD: Handle Undefined Metric Conditions Pre-Aggregation
    # ===========================================================================
    fill_zero_cols = ['Accuracy', 'Sensitivity', 'Contig_Precision', 'BP_Precision', 
                      'F1_Score', 'Mean_Identity', 'Mean_Completeness', 'Contamination_Pct', 
                      'Circularization_Pct', 'TP', 'FP', 'FN', 'Circular_TP']
    
    for col in fill_zero_cols:
        if col in df.columns:
            df[col] = df[col].fillna(0.0)

    # --- Exclude plasmid-free TN rows from metric calculations ----------------
    df_active = df[
        ~((df['Ref_Plasmid_Count'] == 0) & (df['Pred_Plasmid_Count'] == 0))
    ].copy()

    n_excluded = len(df) - len(df_active)
    if n_excluded > 0:
        print(f"Note: {n_excluded} plasmid-free TN row(s) excluded from metric calculations")

    # ===========================================================================
    # 1. Summary by Tool (Macro & Weighted Macro Aggregations)
    # ===========================================================================
    print("\n=== Summary by Tool ===")

    # Helper function for weighted average grouping
    def weighted_mean(group, value_col, weight_col):
        d = group[value_col]
        w = group[weight_col]
        return np.where(w.sum() > 0, np.sum(d * w) / np.sum(w), 0.0).item()

    # Calculate global metrics group by group to allow weighted metrics cleanly
    tool_groups = []
    for tool, group in df_active.groupby('Tool'):
        tool_groups.append({
            'Tool': tool,
            'Samples': len(group),
            'Total_Ref_Plasmids': group['Ref_Plasmid_Count'].sum(),
            'Total_Pred_Plasmids': group['Pred_Plasmid_Count'].sum(),
            'Macro_Accuracy': group['Accuracy'].mean(),
            'Macro_Precision': group['BP_Precision'].mean(),
            'Macro_Sensitivity': group['Sensitivity'].mean(),
            'Macro_F1': group['F1_Score'].mean(),
            # Weighted by Ref_Plasmid_Count to reduce penalty of small/empty artifact samples
            'Weighted_F1': weighted_mean(group, 'F1_Score', 'Ref_Plasmid_Count'),
            'Mean_Circularization': group['Circularization_Pct'].mean(),
            'Mean_Identity': group['Mean_Identity'].mean(),
            'Mean_Completeness': group['Mean_Completeness'].mean(),
            'Mean_Contamination': group['Contamination_Pct'].mean(),
        })
    
    summary = pd.DataFrame(tool_groups).set_index('Tool')
    summary = summary.round(4)
    print(summary.to_string())
    summary.to_csv(os.path.join(ROOT_DIR, 'summary_by_tool.tsv'), sep='\t', index=True)

    # ===========================================================================
    # 2. Summary by Tool × Platform (Macro-Averages)
    # ===========================================================================
    print("\n=== Summary by Tool × Platform ===")

    platform_summary = df_active.groupby(['Tool', 'Platform']).agg(
        Samples              = ('Sample',              'count'),
        Mean_Accuracy        = ('Accuracy',            'mean'),
        Mean_F1              = ('F1_Score',            'mean'),
        Mean_Sensitivity     = ('Sensitivity',         'mean'),
        Mean_BP_Precision    = ('BP_Precision',        'mean'),
        Mean_Contig_Precision= ('Contig_Precision',    'mean'),
        Mean_Completeness    = ('Mean_Completeness',   'mean'),
        Mean_Contamination   = ('Contamination_Pct',   'mean'),
        Mean_Circularization = ('Circularization_Pct', 'mean'),
        Mean_Identity        = ('Mean_Identity',       'mean'),
    ).round(4)

    print(platform_summary.to_string())
    platform_summary.to_csv(os.path.join(ROOT_DIR, 'summary_by_tool_platform.tsv'), sep='\t', index=True)

    # ===========================================================================
    # 3. Summary by Tool × Species
    # ===========================================================================
    print("\n=== Summary by Tool × Species ===")

    species_summary = df_active.groupby(['Tool', 'Species']).agg(
        Samples              = ('Sample',              'count'),
        Mean_Accuracy        = ('Accuracy',            'mean'),
        Mean_F1              = ('F1_Score',            'mean'),
        Mean_Sensitivity     = ('Sensitivity',         'mean'),
        Mean_BP_Precision    = ('BP_Precision',        'mean'),
        Mean_Completeness    = ('Mean_Completeness',   'mean'),
        Mean_Contamination   = ('Contamination_Pct',   'mean'),
        Mean_Circularization = ('Circularization_Pct', 'mean'),
    ).round(4)

    print(species_summary.to_string())
    species_summary.to_csv(os.path.join(ROOT_DIR, 'summary_by_tool_species.tsv'), sep='\t', index=True)

    # ===========================================================================
    # 4. Summary by Tool × Size Category
    # ===========================================================================
    print("\n=== Summary by Size Category ===")

    df_size = df_active[df_active['Size_Category'].notna() & (df_active['Size_Category'] != 'None')].copy()

    if df_size.empty:
        print("WARNING: No rows with valid Size_Category. Check get_ref_size_category.py execution status.")
    else:
        size_order = ['0-1kb', '1-5kb', '5-10kb', '10-100kb', '>100kb']
        df_size['Size_Category'] = pd.Categorical(
            df_size['Size_Category'],
            categories=[s for s in size_order if s in df_size['Size_Category'].unique()],
            ordered=True
        )

        size_summary = df_size.groupby(['Tool', 'Size_Category'], observed=True).agg(
            Samples              = ('Sample',              'count'),
            Ref_Plasmid_Count    = ('Ref_Plasmid_Count',  'sum'),
            Pred_Plasmid_Count   = ('Pred_Plasmid_Count', 'sum'),
            Mean_Accuracy        = ('Accuracy',            'mean'),
            Mean_F1              = ('F1_Score',            'mean'),
            Mean_Sensitivity     = ('Sensitivity',         'mean'),
            Mean_BP_Precision    = ('BP_Precision',        'mean'),
            Mean_Completeness    = ('Mean_Completeness',   'mean'),
            Mean_Contamination   = ('Contamination_Pct',   'mean'),
            Mean_Circularization = ('Circularization_Pct', 'mean'),
        ).round(4)

        print(size_summary.to_string())
        size_summary.to_csv(os.path.join(ROOT_DIR, 'summary_by_size_category.tsv'), sep='\t', index=True)

    # ===========================================================================
    # 5. Per-sample Tool Comparison (Safe pivot table)
    # ===========================================================================
    print("\n=== Per-sample Tool Comparison ===")

    metrics_to_compare = ['Accuracy', 'F1_Score', 'Sensitivity', 'BP_Precision', 
                          'Mean_Completeness', 'Circularization_Pct', 'Contamination_Pct']
    metrics_to_compare = [m for m in metrics_to_compare if m in df_active.columns]

    label_map = {
        'Accuracy': 'Accuracy', 'F1_Score': 'F1', 'Sensitivity': 'Sensitivity',
        'BP_Precision': 'Precision', 'Mean_Completeness': 'Completeness',
        'Circularization_Pct': 'Circularization', 'Contamination_Pct': 'Contamination',
    }

    try:
        comparison = df_active.pivot_table(
            index   = ['Sample', 'Platform', 'Species'],
            columns = 'Tool',
            values  = metrics_to_compare,
            aggfunc = 'mean'
        )

        comparison.columns = [
            f"{label_map.get(metric, metric)}_{tool}"
            for metric, tool in comparison.columns
        ]

        comparison = comparison.fillna(0.0)
        comparison = comparison.reset_index().round(4)
        print(comparison.to_string(index=False))
        comparison.to_csv(os.path.join(ROOT_DIR, 'per_sample_comparison.tsv'), sep='\t', index=False)
    except Exception as e:
        print(f"WARNING: Per-sample comparison failed: {e}")

    # ===========================================================================
    # 6. Quick console summary
    # ===========================================================================
    print("\n=== Files Written ===")
    for fname in ['summary_by_tool.tsv', 'summary_by_tool_platform.tsv', 
                  'summary_by_tool_species.tsv', 'summary_by_size_category.tsv', 'per_sample_comparison.tsv']:
        fpath = os.path.join(ROOT_DIR, fname)
        print(f"  ✓ {fname}" if os.path.isfile(fpath) else f"  ✗ {fname}  (not written)")


if __name__ == '__main__':
    main()