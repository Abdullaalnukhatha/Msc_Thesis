#!/bin/bash
#SBATCH --job-name=assembly_qc
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=24:00:00
#SBATCH --output=/rds/projects/e/elhamsak-mbru-amr/tempdata/assembly_qc_%A_%a.out
#SBATCH --error=/rds/projects/e/elhamsak-mbru-amr/tempdata/assembly_qc_%A_%a.err

set -euo pipefail

# ============================================================
# Assembly QC SLURM array: QUAST + BUSCO + CheckM2
#
# This script DOES NOT create the manifest.
# Create the manifest first, then submit with:
#
#   TOTAL=$(wc -l < /rds/projects/e/elhamsak-mbru-amr/tempdata/assembly_qc_manifest.tsv)
#   sbatch --array=1-${TOTAL}%50 assembly_qc_array.slurm
#
# Each array task processes one assembly path from one manifest line.
# ============================================================

PROJ="/rds/projects/e/elhamsak-mbru-amr"
OUTPUT_DIR="$PROJ/results/assembly_qc"
REFERENCE_DIR="$PROJ/data/ReferenceFASTA"
REFERENCE_MAP="$PROJ/data/ReferenceFASTA/reference_mapping.tsv"
MANIFEST="$PROJ/tempdata/assembly_qc_manifest.tsv"

CHECKM2_DB_FILE="$PROJ/databases/checkm2_db/CheckM2_database/uniref100.KO.1.dmnd"
BUSCO_DOWNLOADS="$PROJ/busco_downloads"
BUSCO_LINEAGE="$BUSCO_DOWNLOADS/lineages/enterobacterales_odb10"
THREADS="${SLURM_CPUS_PER_TASK:-12}"

MINIFORGE="$HOME/miniforge3"
QUAST_ENV="$MINIFORGE/envs/quast_env"
BUSCO_ENV="$MINIFORGE/envs/busco_env"
CHECKM2_ENV="$MINIFORGE/envs/checkm2_env"

mkdir -p "$PROJ/tempdata" "$OUTPUT_DIR"

[[ -s "$MANIFEST" ]] || {
    echo "ERROR: Manifest not found or empty: $MANIFEST" >&2
    exit 1
}

TOTAL="$(wc -l < "$MANIFEST")"
if (( SLURM_ARRAY_TASK_ID < 1 || SLURM_ARRAY_TASK_ID > TOTAL )); then
    echo "ERROR: Array task ID $SLURM_ARRAY_TASK_ID is outside manifest range 1-$TOTAL" >&2
    exit 1
fi

CONTIGS="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")"
[[ -n "$CONTIGS" && -s "$CONTIGS" ]] || {
    echo "ERROR: Assembly missing or empty for task $SLURM_ARRAY_TASK_ID: $CONTIGS" >&2
    exit 1
}

module purge
module load bluebear
module load bear-apps/2024a/live
module load Miniforge3/25.3.0-3

if [[ ! -f "$MINIFORGE/etc/profile.d/conda.sh" ]]; then
    echo "ERROR: Conda initialization script not found: $MINIFORGE/etc/profile.d/conda.sh" >&2
    exit 1
fi
source "$MINIFORGE/etc/profile.d/conda.sh"

for env_path in "$QUAST_ENV" "$BUSCO_ENV" "$CHECKM2_ENV"; do
    [[ -d "$env_path" ]] || {
        echo "ERROR: Conda environment not found: $env_path" >&2
        exit 1
    }
done

if [[ -x "$QUAST_ENV/bin/quast.py" ]]; then
    QUAST_COMMAND="quast.py"
elif [[ -x "$QUAST_ENV/bin/quast" ]]; then
    QUAST_COMMAND="quast"
else
    echo "ERROR: QUAST not found in: $QUAST_ENV" >&2
    exit 1
fi

[[ -x "$BUSCO_ENV/bin/python" ]] || {
    echo "ERROR: BUSCO Python not found: $BUSCO_ENV/bin/python" >&2
    exit 1
}
[[ -f "$BUSCO_ENV/bin/busco" ]] || {
    echo "ERROR: BUSCO launcher not found: $BUSCO_ENV/bin/busco" >&2
    exit 1
}
[[ -x "$BUSCO_ENV/bin/hmmsearch" ]] || {
    echo "ERROR: HMMER hmmsearch not found: $BUSCO_ENV/bin/hmmsearch" >&2
    exit 1
}
"$BUSCO_ENV/bin/python" -c "import busco" || {
    echo "ERROR: BUSCO Python module cannot be imported from: $BUSCO_ENV" >&2
    exit 1
}
[[ -x "$CHECKM2_ENV/bin/checkm2" ]] || {
    echo "ERROR: CheckM2 not found: $CHECKM2_ENV/bin/checkm2" >&2
    exit 1
}

[[ -f "$REFERENCE_MAP" ]] || {
    echo "ERROR: Reference mapping not found: $REFERENCE_MAP" >&2
    exit 1
}
[[ -d "$REFERENCE_DIR" ]] || {
    echo "ERROR: Reference directory not found: $REFERENCE_DIR" >&2
    exit 1
}
[[ -d "$BUSCO_LINEAGE" ]] || {
    echo "ERROR: BUSCO lineage not found: $BUSCO_LINEAGE" >&2
    exit 1
}
[[ -s "$CHECKM2_DB_FILE" ]] || {
    echo "ERROR: CheckM2 database file not found: $CHECKM2_DB_FILE" >&2
    exit 1
}

FILENAME="$(basename "$CONTIGS")"
if [[ "$FILENAME" == "assembly.fasta" ]]; then
    SAMPLE="$(basename "$(dirname "$CONTIGS")")"
    SPECIES="$(basename "$(dirname "$(dirname "$CONTIGS")")")"
    PLATFORM="$(basename "$(dirname "$(dirname "$(dirname "$CONTIGS")")")")"
else
    SAMPLE="${FILENAME%_contigs.fasta}"
    SPECIES="$(basename "$(dirname "$CONTIGS")")"
    PLATFORM="$(basename "$(dirname "$(dirname "$CONTIGS")")")"
fi

SAMPLE_OUT="$OUTPUT_DIR/$PLATFORM/$SPECIES/$SAMPLE"

# Every submitted task regenerates this sample's QC outputs.
rm -rf "$SAMPLE_OUT"
mkdir -p "$SAMPLE_OUT"

 echo "=================================================================="
 echo " Assembly QC array task"
 echo "=================================================================="
 echo "Array job/task : ${SLURM_ARRAY_JOB_ID:-unknown}/${SLURM_ARRAY_TASK_ID}"
 echo "Manifest line  : $SLURM_ARRAY_TASK_ID of $TOTAL"
 echo "Sample         : $SAMPLE"
 echo "Species        : $SPECIES"
 echo "Platform       : $PLATFORM"
 echo "Assembly       : $CONTIGS"
 echo "Output         : $SAMPLE_OUT"
 echo "Threads        : $THREADS"
 echo "Started        : $(date)"
 echo "=================================================================="

# ----------------------------- QUAST -----------------------------
REF_ID="$(awk -F '\t' -v sample="$SAMPLE" '$1 == sample {print $2; exit}' "$REFERENCE_MAP")"
REF_FASTA=""

if [[ -n "$REF_ID" ]]; then
    # Prefer a non-empty whole-genome reference. Do not use plasmid-only FASTA files.
    REF_FASTA="$(find "$REFERENCE_DIR" -type f -size +0c \
        -name "${REF_ID}*_genomic.fna" \
        ! -name "*_plasmid.fna" \
        -print -quit)"

    # Fallback for valid non-empty genomic FASTA files with a different extension.
    if [[ -z "$REF_FASTA" ]]; then
        REF_FASTA="$(find "$REFERENCE_DIR" -type f -size +0c \
            \( -name "${REF_ID}*.fasta" -o -name "${REF_ID}*.fa" \) \
            ! -name "*_plasmid.fasta" ! -name "*_plasmid.fa" \
            -print -quit)"
    fi
fi

if [[ -n "$REF_FASTA" && -s "$REF_FASTA" ]]; then
    echo "QUAST reference: $REF_FASTA"
    conda run --no-capture-output -p "$QUAST_ENV" "$QUAST_COMMAND" "$CONTIGS" \
        -r "$REF_FASTA" \
        --output-dir "$SAMPLE_OUT/quast" \
        --threads "$THREADS" \
        --min-contig 500
else
    echo "WARNING: No reference found for $SAMPLE; running QUAST without one."
    conda run --no-capture-output -p "$QUAST_ENV" "$QUAST_COMMAND" "$CONTIGS" \
        --output-dir "$SAMPLE_OUT/quast" \
        --threads "$THREADS" \
        --min-contig 500
fi

# ----------------------------- BUSCO -----------------------------
PATH="$BUSCO_ENV/bin:$PATH" \
"$BUSCO_ENV/bin/python" "$BUSCO_ENV/bin/busco" \
    -i "$CONTIGS" \
    -o busco \
    --out_path "$SAMPLE_OUT" \
    -l "$BUSCO_LINEAGE" \
    -m genome \
    --cpu "$THREADS" \
    --offline \
    -f

# ---------------------------- CheckM2 ----------------------------
conda run --no-capture-output -p "$CHECKM2_ENV" checkm2 predict \
    --input "$CONTIGS" \
    --output-directory "$SAMPLE_OUT/checkm2" \
    --database_path "$CHECKM2_DB_FILE" \
    --threads "$THREADS" \
    --force

echo "=================================================================="
echo "Completed: $SAMPLE"
echo "Finished : $(date)"
echo "=================================================================="
