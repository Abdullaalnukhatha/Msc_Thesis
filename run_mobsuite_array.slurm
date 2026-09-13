#!/bin/bash
#SBATCH --job-name=mobsuite_array
#SBATCH --partition=emerald-shared
#SBATCH --time=08:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-120%20
#SBATCH --output=/rds/projects/e/elhamsak-mbru-amr/tempdata/mobsuite_%A_%a.out

set -euo pipefail

module purge
module load bluebear
module load bear-apps/2024a/live

source /rds/homes/a/axa2273/miniforge3/etc/profile.d/conda.sh
conda activate /rds/projects/e/elhamsak-mbru-amr/env/mobsuite_env

MANIFEST="/rds/projects/e/elhamsak-mbru-amr/results/Best_Assemblies/master_manifest.tsv"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

if [ ! -f "$MANIFEST" ]; then
    echo "ERROR: Manifest not found: $MANIFEST"
    exit 1
fi

# get row safely
ROW=$(tail -n +2 "$MANIFEST" | sed -n "${SLURM_ARRAY_TASK_ID}p" || true)
ROW=$(echo "$ROW" | tr -d '\r')

if [ -z "$ROW" ]; then
    echo "ERROR: No manifest row found for task $SLURM_ARRAY_TASK_ID"
    exit 1
fi

# master_manifest has 6 columns, but we only want the first 5
SAMPLE=$(echo "$ROW"   | cut -f1)
PLATFORM=$(echo "$ROW" | cut -f2)
SPECIES=$(echo "$ROW"  | cut -f3)
GROUP=$(echo "$ROW"    | cut -f4)
FASTA=$(echo "$ROW"    | cut -f5)

OUTDIR="/rds/projects/e/elhamsak-mbru-amr/results/plasmid_benchmark/MOB-suite/$PLATFORM/$SPECIES/$SAMPLE"

echo "Parsed values:"
echo "  SAMPLE   = [$SAMPLE]"
echo "  PLATFORM = [$PLATFORM]"
echo "  SPECIES  = [$SPECIES]"
echo "  GROUP    = [$GROUP]"
echo "  FASTA    = [$FASTA]"
echo "  OUTDIR   = [$OUTDIR]"
echo

if [ ! -f "$FASTA" ]; then
    echo "ERROR: FASTA not found: $FASTA"
    exit 1
fi

echo "Processing $SAMPLE"

# MOB-suite wants OUTDIR itself to not exist
rm -rf "$OUTDIR"
mkdir -p "$(dirname "$OUTDIR")"

mob_recon \
  --infile "$FASTA" \
  --outdir "$OUTDIR" \
  --sample_id "$SAMPLE" \
  --num_threads "$THREADS" \
  --run_overhang

cat "$OUTDIR"/plasmid_*.fasta 2>/dev/null > "$OUTDIR/${SAMPLE}_plasmids.fasta" || touch "$OUTDIR/${SAMPLE}_plasmids.fasta"

echo "Plasmids found: $(grep -c '^>' "$OUTDIR/${SAMPLE}_plasmids.fasta" 2>/dev/null || echo 0)"
echo "DONE: $SAMPLE"
