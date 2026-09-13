#!/usr/bin/env python3
# select_best_assemblies.py
# Select best assemblies based on N50 + CheckM2 completeness

import os
import sys
import shutil
import argparse
import csv
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Select best assemblies based on N50 + CheckM2 completeness"
    )
    project_root = Path("/rds/projects/e/elhamsak-mbru-amr")

    parser.add_argument(
        "--qc-dir",
        default=str(project_root / "results" / "assembly_qc"),
        help="QC results directory",
    )
    parser.add_argument(
        "--assembly-dir",
        default=str(project_root / "results" / "assemblies"),
        help="Assembly directory",
    )
    parser.add_argument(
        "--output-dir",
        default=str(project_root / "results" / "Best_Assemblies"),
        help="Output directory",
    )
    parser.add_argument("--n",            type=int, default=60,           help="Number of assemblies to select per group")
    parser.add_argument("--w-completeness", type=float, default=0.7,      help="Weight for completeness score")
    parser.add_argument("--w-n50",          type=float, default=0.3,      help="Weight for N50 score")
    parser.add_argument("--contamination-threshold", type=float, default=5.0,
                        help="Maximum allowed contamination (%)")
    return parser.parse_args()


def parse_quast_report(report_path):
    """
    Extract N50 from QUAST report.tsv.
    """
    with open(report_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if parts[0] == "N50":
                return float(parts[1])
    return None


def parse_checkm2_report(report_path):
    """
    Extract completeness and contamination from CheckM2 quality_report.tsv.
    """
    with open(report_path) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            return float(row['Completeness']), float(row['Contamination'])
    return None, None


def collect_metrics(qc_dir, platform, species):
    """
    Collect QC metrics for all samples in a platform/species group.
    """
    species_dir = Path(qc_dir) / platform / species
    samples = []

    for sample_dir in sorted(species_dir.iterdir()):
        if not sample_dir.is_dir():
            continue
        sample = sample_dir.name

        quast_report  = sample_dir / "quast" / "report.tsv"
        checkm2_report = sample_dir / "checkm2" / "quality_report.tsv"

        if not quast_report.exists():
            print(f"  WARNING: QUAST report missing for {sample}, skipping")
            continue
        if not checkm2_report.exists():
            print(f"  WARNING: CheckM2 report missing for {sample}, skipping")
            continue

        n50 = parse_quast_report(quast_report)
        completeness, contamination = parse_checkm2_report(checkm2_report)

        if n50 is None or completeness is None or contamination is None:
            print(f"  WARNING: Missing metrics for {sample}, skipping")
            continue

        samples.append({
            "sample":        sample,
            "n50":           n50,
            "completeness":  completeness,
            "contamination": contamination,
        })

    return samples


def filter_and_score(samples, w_completeness, w_n50, contamination_threshold):
    """
    Filter by contamination and compute normalized weighted score.
    """

    # Filter contamination
    filtered = [s for s in samples if s["contamination"] <= contamination_threshold]
    n_removed = len(samples) - len(filtered)
    if n_removed > 0:
        print(f"  Removed {n_removed} samples with contamination > {contamination_threshold}%")

    if not filtered:
        return []

    # Normalize N50 and completeness to [0, 1]
    max_n50  = max(s["n50"]          for s in filtered)
    max_comp = max(s["completeness"] for s in filtered)

    for s in filtered:
        norm_n50  = s["n50"]          / max_n50  if max_n50  > 0 else 0
        norm_comp = s["completeness"] / max_comp if max_comp > 0 else 0
        s["score"] = w_completeness * norm_comp + w_n50 * norm_n50

    # Sort by score descending
    filtered.sort(key=lambda x: x["score"], reverse=True)
    return filtered


def write_summary(selected, output_dir, platform, species):
    """
    Write selected samples summary TSV.
    """
    summary_path = Path(output_dir) / f"{platform}_{species}_selected.tsv"
    with open(summary_path, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=[
            "SAMPLE", "N50", "COMPLETENESS", "CONTAMINATION", "SCORE"
        ], delimiter='\t')
        writer.writeheader()
        for s in selected:
            writer.writerow({
                "SAMPLE":        s["sample"],
                "N50":           f"{s['n50']:.0f}",
                "COMPLETENESS":  f"{s['completeness']:.2f}",
                "CONTAMINATION": f"{s['contamination']:.2f}",
                "SCORE":         f"{s['score']:.4f}",
            })
    return summary_path


def copy_assemblies(selected, assembly_dir, output_dir, platform, species):
    """
    Copy selected assembly FASTA files and any matching assembly_info files
    to the output directory.
    """
    src_dir = Path(assembly_dir) / platform / species
    dest_dir = Path(output_dir) / platform / species
    dest_dir.mkdir(parents=True, exist_ok=True)

    missing = []

    for s in selected:
        sample = s["sample"]

        # Copy contigs FASTA
        fasta_files = list(src_dir.glob(f"{sample}_contigs.fasta"))
        if fasta_files:
            shutil.copy2(fasta_files[0], dest_dir / fasta_files[0].name)
        else:
            print(f"  WARNING: Assembly not found for {sample}")
            missing.append(sample)

        # Copy all matching assembly_info files
        info_files = list(src_dir.glob("*assembly_info.txt"))
        for info_file in info_files:
            shutil.copy2(info_file, dest_dir / info_file.name)

    return dest_dir, missing

def main():
    args = parse_args()

    print("=== Best Assembly Selection ===")
    print(f"QC directory:             {args.qc_dir}")
    print(f"Assembly directory:       {args.assembly_dir}")
    print(f"Output directory:         {args.output_dir}")
    print(f"Default N per group:      {args.n}")
    print(f"Weights:                  completeness={args.w_completeness}, N50={args.w_n50}")
    print(f"Contamination threshold:  >{args.contamination_threshold}%")
    print()

    os.makedirs(args.output_dir, exist_ok=True)

    qc_path = Path(args.qc_dir)
    if not qc_path.exists():
        print(f"ERROR: QC directory not found: {args.qc_dir}")
        sys.exit(1)

    summary_results = []

    platforms = sorted([p.name for p in qc_path.iterdir() if p.is_dir()])
    if not platforms:
        print("ERROR: No platform directories found in QC directory")
        sys.exit(1)

    for platform in platforms:
        platform_dir = qc_path / platform
        species_list = sorted([s.name for s in platform_dir.iterdir() if s.is_dir()])

        for species in species_list:
            print(f"=== Processing: {platform} / {species} ===")

            # Collect metrics
            samples = collect_metrics(args.qc_dir, platform, species)
            print(f"  Total samples found:    {len(samples)}")

            if not samples:
                print("  WARNING: No samples with complete metrics, skipping")
                print()
                continue

            # Filter and score
            scored = filter_and_score(
                samples,
                args.w_completeness,
                args.w_n50,
                args.contamination_threshold
            )
            print(f"  After filtering:        {len(scored)} samples remaining")

            if not scored:
                print("  WARNING: No samples passed filters, skipping")
                print()
                continue

            # Select top N
            # For groups with fewer samples than args.n, select all that pass filters
            n_select = min(args.n, len(scored))
            selected = scored[:n_select]
            print(f"  Selecting top:          {n_select} samples")

            # Write summary TSV
            summary_path = write_summary(selected, args.output_dir, platform, species)
            print(f"  Summary written:        {summary_path}")

            # Copy assemblies
            dest_dir, missing = copy_assemblies(
                selected, args.assembly_dir, args.output_dir, platform, species
            )
            print(f"  Assemblies copied to:   {dest_dir}")
            if missing:
                print(f"  WARNING: {len(missing)} assemblies not found: {', '.join(missing)}")

            summary_results.append((platform, species, n_select, summary_path))
            print()

    # Final summary
    print("=== Selection Complete ===")
    print()
    print("Summary:")
    for platform, species, n, path in summary_results:
        print(f"  {platform}/{species}: {n} samples selected → {path}")


if __name__ == "__main__":
    main()
