#!/bin/bash
#SBATCH --job-name=PacBio_filter
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/rds/projects/e/elhamsak-mbru-amr/tempdata/PacBio_filter_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-mbru-amr/tempdata/PacBio_filter_%j.err

set -euo pipefail

# ------------------------------------------------------------------
# Modules and Conda environment
# ------------------------------------------------------------------

module purge
module load bluebear
module load bear-apps/2024a/live
module load Miniforge3/25.3.0-3

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /rds/homes/a/axa2273/miniforge3/envs/porechop_env

# ------------------------------------------------------------------
# Configuration
#
# Expected input structure:
#
# data/FASTQ/
# ├── ecoli/
# │   ├── ONT/
# │   └── PacBio/
# └── kpneumoniae/
#     ├── ONT/
#     └── PacBio/
#
# Output structure:
#
# results/ProcessedFASTQ/PacBio/
# ├── ecoli/
# └── kpneumoniae/
# ------------------------------------------------------------------

PROJ="/rds/projects/e/elhamsak-mbru-amr"

INPUT_ROOT="${1:-$PROJ/data/FASTQ}"
OUTPUT_ROOT="${2:-$PROJ/results/ProcessedFASTQ/PacBio}"

MIN_LENGTH="${3:-500}"
MIN_QUALITY="${4:-7}"
THREADS="${SLURM_CPUS_PER_TASK:-12}"

SPECIES_LIST=("ecoli" "kpneumoniae")

mkdir -p "$OUTPUT_ROOT"
mkdir -p "$PROJ/tempdata"

# ------------------------------------------------------------------
# Check required commands
# ------------------------------------------------------------------

for cmd in seqkit gzip bc find; do
    if ! command -v "$cmd" >/dev/null 2>&1; then
        echo "ERROR: Required command not found: $cmd" >&2
        exit 1
    fi
done

# ------------------------------------------------------------------
# Temporary-file cleanup
# ------------------------------------------------------------------

CURRENT_TMP=""

cleanup() {
    if [[ -n "$CURRENT_TMP" && -f "$CURRENT_TMP" ]]; then
        echo "Removing unfinished temporary file: $CURRENT_TMP"
        rm -f "$CURRENT_TMP"
    fi
}

trap cleanup EXIT INT TERM

# ------------------------------------------------------------------
# Job information
# ------------------------------------------------------------------

echo "=================================================================="
echo " PacBio read filtering"
echo "=================================================================="
echo "Job ID       : ${SLURM_JOB_ID:-N/A}"
echo "Node         : $(hostname)"
echo "Input root   : $INPUT_ROOT"
echo "Output root  : $OUTPUT_ROOT"
echo "Min length   : $MIN_LENGTH bp"
echo "Min quality  : Q$MIN_QUALITY"
echo "Threads      : $THREADS"
echo "Started      : $(date)"
echo "=================================================================="

if [[ ! -d "$INPUT_ROOT" ]]; then
    echo "ERROR: FASTQ input root does not exist:"
    echo "$INPUT_ROOT"
    exit 1
fi

# ------------------------------------------------------------------
# Process each species
# ------------------------------------------------------------------

for SPECIES in "${SPECIES_LIST[@]}"; do

    SPECIES_INPUT="$INPUT_ROOT/$SPECIES/PacBio"
    SPECIES_OUTPUT="$OUTPUT_ROOT/$SPECIES"

    echo
    echo "=================================================================="
    echo "Processing species: $SPECIES"
    echo "Input             : $SPECIES_INPUT"
    echo "Output            : $SPECIES_OUTPUT"
    echo "=================================================================="

    if [[ ! -d "$SPECIES_INPUT" ]]; then
        echo "WARNING: PacBio directory not found, skipping:"
        echo "$SPECIES_INPUT"
        continue
    fi

    mkdir -p "$SPECIES_OUTPUT"

    mapfile -t FASTQ_FILES < <(
        find "$SPECIES_INPUT" -type f \
            \( -name "*.fastq.gz" -o -name "*.fq.gz" \) \
            | sort
    )

    TOTAL_FILES=${#FASTQ_FILES[@]}

    if (( TOTAL_FILES == 0 )); then
        echo "WARNING: No compressed FASTQ files found under:"
        echo "$SPECIES_INPUT"
        continue
    fi

    echo "Found $TOTAL_FILES PacBio FASTQ files."

    STATS_FILE="$SPECIES_OUTPUT/${SPECIES}_processing_stats.tsv"

    printf \
        "Sample\tInput_Reads\tInput_Bases\tFiltered_Reads\tFiltered_Bases\tRetention_Rate\n" \
        > "$STATS_FILE"

    CURRENT=0
    COMPLETED=0
    FAILED=0

    for FASTQ in "${FASTQ_FILES[@]}"; do

        ((CURRENT += 1))

        FILENAME=$(basename "$FASTQ")
        SAMPLE=$(echo "$FILENAME" | sed -E 's/\.(fastq|fq)\.gz$//')

        OUT_FASTQ="$SPECIES_OUTPUT/${SAMPLE}_processed.fastq.gz"
        TMP_FASTQ="$SPECIES_OUTPUT/.${SAMPLE}_processed.tmp.fastq.gz"

        CURRENT_TMP="$TMP_FASTQ"

        echo
        echo "[$CURRENT/$TOTAL_FILES] Processing: $SAMPLE"
        echo "Input : $FASTQ"
        echo "Output: $OUT_FASTQ"

        rm -f "$TMP_FASTQ"

        # Validate compressed input.
        if ! gzip -t "$FASTQ" 2>/dev/null; then
            echo "WARNING: Corrupted gzip input; skipping: $FASTQ"
            ((FAILED += 1))
            CURRENT_TMP=""
            continue
        fi

        # Get input statistics.
        INPUT_STATS=$(seqkit stats -T "$FASTQ" | tail -n 1)

        INPUT_READS=$(echo "$INPUT_STATS" | cut -f4 | tr -d ',')
        INPUT_BASES=$(echo "$INPUT_STATS" | cut -f5 | tr -d ',')

        if [[ ! "$INPUT_READS" =~ ^[0-9]+$ ]]; then
            echo "WARNING: Could not determine input read count for $SAMPLE"
            ((FAILED += 1))
            CURRENT_TMP=""
            continue
        fi

        # Filter by read length and average quality.
        if ! seqkit seq \
            -m "$MIN_LENGTH" \
            -Q "$MIN_QUALITY" \
            -j "$THREADS" \
            "$FASTQ" \
            -o "$TMP_FASTQ"; then

            echo "WARNING: SeqKit filtering failed for $SAMPLE"
            rm -f "$TMP_FASTQ"
            ((FAILED += 1))
            CURRENT_TMP=""
            continue
        fi

        # Validate output before renaming.
        if ! gzip -t "$TMP_FASTQ" 2>/dev/null; then
            echo "WARNING: Filtered output is corrupted for $SAMPLE"
            rm -f "$TMP_FASTQ"
            ((FAILED += 1))
            CURRENT_TMP=""
            continue
        fi

        mv -f "$TMP_FASTQ" "$OUT_FASTQ"
        CURRENT_TMP=""

        # Get output statistics.
        OUTPUT_STATS=$(seqkit stats -T "$OUT_FASTQ" | tail -n 1)

        OUTPUT_READS=$(echo "$OUTPUT_STATS" | cut -f4 | tr -d ',')
        OUTPUT_BASES=$(echo "$OUTPUT_STATS" | cut -f5 | tr -d ',')

        if [[ ! "$OUTPUT_READS" =~ ^[0-9]+$ ]]; then
            echo "WARNING: Could not determine output read count for $SAMPLE"
            ((FAILED += 1))
            continue
        fi

        if (( INPUT_READS > 0 )); then
            RETENTION=$(echo "scale=2; $OUTPUT_READS * 100 / $INPUT_READS" | bc)
        else
            RETENTION="0.00"
        fi

        printf "%s\t%s\t%s\t%s\t%s\t%s%%\n" \
            "$SAMPLE" \
            "$INPUT_READS" \
            "$INPUT_BASES" \
            "$OUTPUT_READS" \
            "$OUTPUT_BASES" \
            "$RETENTION" \
            >> "$STATS_FILE"

        echo "Retained: $OUTPUT_READS / $INPUT_READS reads (${RETENTION}%)"

        ((COMPLETED += 1))
    done

    echo
    echo "------------------------------------------------------------------"
    echo "Species completed : $SPECIES"
    echo "Total inputs      : $TOTAL_FILES"
    echo "Completed         : $COMPLETED"
    echo "Failed/skipped    : $FAILED"
    echo "Statistics        : $STATS_FILE"
    echo "------------------------------------------------------------------"

done

echo
echo "=================================================================="
echo "PacBio filtering completed: $(date)"
echo "=================================================================="
