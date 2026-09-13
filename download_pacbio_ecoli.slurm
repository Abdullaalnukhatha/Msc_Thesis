#!/bin/bash
#SBATCH --job-name=Ecoli_Pacbio_download
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --output=/rds/projects/e/elhamsak-mbru-amr/tempdata/download_%j_ecoli_pacbio.out

set -euo pipefail

module purge
module load bluebear
module load bear-apps/2024a/live
module load Miniforge3/25.3.0-3

source /rds/homes/a/axa2273/miniforge3/etc/profile.d/conda.sh
conda activate sra

export PATH="$HOME/.conda/envs/sra/bin:$PATH"
echo "prefetch:     $(which prefetch || echo NOT_FOUND)"
echo "fasterq-dump: $(which fasterq-dump)"

# --- Project base (avoid home quota) ---
PROJ="/rds/projects/e/elhamsak-mbru-amr"
mkdir -p "$PROJ/tempdata" "$PROJ/scratch"

# Job scratch in PROJECT space (NOT $HOME)
SCRATCH_BASE="$PROJ/scratch"
export TMPDIR="$SCRATCH_BASE/sra_tmp/$SLURM_JOB_ID"
mkdir -p "$TMPDIR"

# Put SRA downloads/cache in PROJECT scratch
SRA_CACHE="$SCRATCH_BASE/sra_cache"
SRA_STORE="$SCRATCH_BASE/sra_store/$SLURM_JOB_ID"
mkdir -p "$SRA_CACHE" "$SRA_STORE"

# SRA toolkit config location (keep in job TMPDIR)
export VDB_CONFIG="$TMPDIR/vdb-config"
export NCBI_SETTINGS="$VDB_CONFIG"
mkdir -p "$VDB_CONFIG"

# (Optional) cache hint
export VDB_CACHE_ROOT="$SRA_CACHE"

LIST="best_pacbio_ecoli_srr_list.txt"

# Output FASTQs into PROJECT space (recommended)
OUTDIR="$PROJ/data/FASTQ/ecoli/PacBio"
mkdir -p "$OUTDIR"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# Clean list (remove empty + comment lines)
CLEAN_LIST="$(mktemp)"
grep -v '^\s*$' "$LIST" | grep -v '^\s*#' > "$CLEAN_LIST"

TOTAL=$(wc -l < "$CLEAN_LIST" | tr -d ' ')
i=0

download_one () {
  local srr="$1"

  # Skip if already downloaded
  if ls "${OUTDIR}/${srr}"*.fastq.gz >/dev/null 2>&1; then
    echo "[SKIP] $srr already exists"
    return 0
  fi

  echo "========================================="
  echo "Downloading: $srr"
  echo "========================================="

  for attempt in 1 2 3; do
    echo "[INFO] $srr attempt $attempt/3"

    # 1) Download .sra locally using prefetch
    if prefetch "$srr" --max-size 0 --output-directory "$SRA_STORE" 2>&1 | tee -a "${OUTDIR}/${srr}_download.log"; then
      :
    else
      echo "[WARN] prefetch failed for $srr on attempt $attempt" | tee -a "${OUTDIR}/${srr}_download.log"
      sleep 60
      continue
    fi

    # Locate downloaded .sra file
    local sra_path="${SRA_STORE}/${srr}/${srr}.sra"
    if [ ! -f "$sra_path" ]; then
      sra_path="$(find "$SRA_STORE" -maxdepth 3 -name "${srr}.sra" 2>/dev/null | head -n 1 || true)"
    fi

    if [ -z "${sra_path:-}" ] || [ ! -f "${sra_path:-}" ]; then
      echo "[WARN] .sra not found after prefetch for $srr (attempt $attempt)" | tee -a "${OUTDIR}/${srr}_download.log"
      sleep 60
      continue
    fi

    # 2) Convert locally -> fastq
    if fasterq-dump "$sra_path" \
        --outdir "$OUTDIR" \
        --split-spot \
        --skip-technical \
        -e "$THREADS" \
        2>&1 | tee -a "${OUTDIR}/${srr}_download.log"
    then
      # gzip only if FASTQ exists
      if ls "${OUTDIR}/${srr}"*.fastq >/dev/null 2>&1; then
        gzip -f "${OUTDIR}/${srr}"*.fastq
      fi
      echo "[OK] $srr"

      # remove .sra to save space
      rm -f "$sra_path" || true
      return 0
    fi

    echo "[WARN] fasterq-dump failed for $srr on attempt $attempt" | tee -a "${OUTDIR}/${srr}_download.log"
    sleep 60
  done

  echo "[FAILED] $srr after 3 attempts (skipping)" | tee -a "${OUTDIR}/${srr}_download.log"
  return 0
}

while read -r srr; do
  i=$((i+1))
  echo ""
  echo "Progress: [$i/$TOTAL]"
  download_one "$srr"
done < "$CLEAN_LIST"

rm -f "$CLEAN_LIST"

# Cleanup job-specific SRA store (keeps cache for performance)
rm -rf "$SRA_STORE" || true

echo ""
echo "========================================="
echo "DONE. Output directory: $OUTDIR"
echo "Tempdata (logs): $PROJ/tempdata"
echo "Scratch base: $SCRATCH_BASE"
echo "========================================="
