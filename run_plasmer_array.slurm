#!/bin/bash
#SBATCH --job-name=plasmer_array
#SBATCH --partition=emerald-shared
#SBATCH --time=08:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=4
#SBATCH --array=1-120%20
#SBATCH --output=/rds/projects/e/elhamsak-mbru-amr/tempdata/plasmer_%A_%a.out
#SBATCH --error=/rds/projects/e/elhamsak-mbru-amr/tempdata/plasmer_%A_%a.err

set -euo pipefail
export LC_ALL=C

PROJ="/rds/projects/e/elhamsak-mbru-amr"
MANIFEST="$PROJ/results/Best_Assemblies/master_manifest.tsv"
OUTPUT_DIR="$PROJ/results/plasmid_benchmark/Plasmer"
DB="$PROJ/databases/PlasmerDB"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

module purge
module load bluebear
module load bear-apps/2024a/live

source /rds/homes/a/axa2273/miniforge3/etc/profile.d/conda.sh
conda activate "$PROJ/env/plasmer_env"

echo "========================================="
echo "Plasmer array job"
echo "Job ID: ${SLURM_JOB_ID:-NA}"
echo "Array job ID: ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID:-NA}"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "Manifest: $MANIFEST"
echo "Output root: $OUTPUT_DIR"
echo "Database: $DB"
echo "Threads: $THREADS"
echo "========================================="

# ------------------------------------------------------------
# Validate SLURM array task
# ------------------------------------------------------------

if [ -z "${SLURM_ARRAY_TASK_ID:-}" ]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not set."
    echo "Submit this script using sbatch, not bash."
    exit 1
fi

# ------------------------------------------------------------
# Validate required files, directories, and command
# ------------------------------------------------------------

if [ ! -f "$MANIFEST" ]; then
    echo "ERROR: Manifest not found:"
    echo "$MANIFEST"
    exit 1
fi

if [ ! -d "$DB" ]; then
    echo "ERROR: Plasmer database directory not found:"
    echo "$DB"
    exit 1
fi

if ! command -v Plasmer >/dev/null 2>&1; then
    echo "ERROR: Plasmer command not found in environment:"
    echo "$PROJ/env/plasmer_env"
    exit 1
fi

# ------------------------------------------------------------
# Read one manifest row
# Header is skipped with tail -n +2
# ------------------------------------------------------------

ROW=$(
    tail -n +2 "$MANIFEST" |
    sed -n "${SLURM_ARRAY_TASK_ID}p" |
    tr -d '\r'
)

if [ -z "$ROW" ]; then
    echo "ERROR: No manifest row found for task ${SLURM_ARRAY_TASK_ID}."
    echo "The array range may be larger than the number of manifest records."
    exit 1
fi

IFS=$'\t' read -r \
    SAMPLE \
    PLATFORM \
    SPECIES \
    GROUP \
    FASTA \
    SOURCE_QC_TSV <<< "$ROW"

# ------------------------------------------------------------
# Validate parsed manifest fields
# ------------------------------------------------------------

if [ -z "${SAMPLE:-}" ] || \
   [ -z "${PLATFORM:-}" ] || \
   [ -z "${SPECIES:-}" ] || \
   [ -z "${FASTA:-}" ]; then
    echo "ERROR: Required manifest fields are missing."
    echo "Manifest row:"
    printf '%s\n' "$ROW"
    exit 1
fi

SAMPLE_OUT="$OUTPUT_DIR/$PLATFORM/$SPECIES/$SAMPLE"
FINAL_FASTA="$SAMPLE_OUT/${SAMPLE}_plasmids.fasta"
PRED_FASTA="$SAMPLE_OUT/results/${SAMPLE}.plasmer.predPlasmids.fa"
STATUS_FILE="$SAMPLE_OUT/plasmer_status.tsv"

echo ""
echo "Parsed values:"
echo "  SAMPLE        = [$SAMPLE]"
echo "  PLATFORM      = [$PLATFORM]"
echo "  SPECIES       = [$SPECIES]"
echo "  GROUP         = [$GROUP]"
echo "  FASTA         = [$FASTA]"
echo "  SOURCE_QC_TSV = [$SOURCE_QC_TSV]"
echo "  OUTDIR        = [$SAMPLE_OUT]"
echo ""

if [ ! -s "$FASTA" ]; then
    echo "ERROR: Input FASTA not found or empty:"
    echo "$FASTA"
    exit 1
fi

# ------------------------------------------------------------
# Skip samples that already completed
#
# The second status-file row must contain:
#   Status = OK or NO_PRED_FASTA
#   Exit_Code = 0
#
# NO_PRED_FASTA is considered complete because Plasmer finished
# successfully but did not create a predicted plasmid FASTA.
# ------------------------------------------------------------

if [ -s "$STATUS_FILE" ]; then
    EXISTING_STATUS=$(
        awk -F'\t' '
            NR == 2 {
                print $4
                exit
            }
        ' "$STATUS_FILE"
    )

    EXISTING_EXIT_CODE=$(
        awk -F'\t' '
            NR == 2 {
                print $5
                exit
            }
        ' "$STATUS_FILE"
    )

    if { [ "$EXISTING_STATUS" = "OK" ] || \
         [ "$EXISTING_STATUS" = "NO_PRED_FASTA" ]; } && \
       [ "$EXISTING_EXIT_CODE" = "0" ]; then

        echo "SKIPPING: Existing completed Plasmer result found."
        echo "Sample: $SAMPLE"
        echo "Existing status: $EXISTING_STATUS"
        echo "Status file: $STATUS_FILE"
        echo "Final FASTA: $FINAL_FASTA"
        echo "Finished: $(date)"
        exit 0
    fi

    echo "Existing result is incomplete or failed."
    echo "Existing status: ${EXISTING_STATUS:-UNKNOWN}"
    echo "Existing exit code: ${EXISTING_EXIT_CODE:-UNKNOWN}"
    echo "The existing sample output will be deleted and rerun."
fi

# ------------------------------------------------------------
# Remove failed or incomplete output before rerunning
# ------------------------------------------------------------

rm -rf "$SAMPLE_OUT"
mkdir -p "$SAMPLE_OUT"

echo ""
echo "Running Plasmer..."
echo ""

set +e

Plasmer \
    -g "$FASTA" \
    -p "$SAMPLE" \
    -o "$SAMPLE_OUT" \
    -d "$DB" \
    -t "$THREADS"

PLASMER_EXIT=$?

set -e

echo ""
echo "Checking Plasmer output..."
echo ""

# ------------------------------------------------------------
# Record a failed Plasmer run
# ------------------------------------------------------------

if [ "$PLASMER_EXIT" -ne 0 ]; then
    echo "WARNING: Plasmer failed for $SAMPLE."
    echo "Exit code: $PLASMER_EXIT"
    echo "Recording the sample as failed with 0 plasmids."

    touch "$FINAL_FASTA"

    {
        echo -e "Sample\tPlatform\tSpecies\tStatus\tExit_Code\tFASTA\tPlasmid_Count\tFinal_FASTA"
        echo -e "$SAMPLE\t$PLATFORM\t$SPECIES\tFAILED\t$PLASMER_EXIT\t$FASTA\t0\t$FINAL_FASTA"
    } > "$STATUS_FILE"

    echo ""
    echo "Plasmids found: 0"
    echo "Final FASTA: $FINAL_FASTA"
    echo "Status file: $STATUS_FILE"
    echo "DONE_WITH_FAILURE_RECORDED: $SAMPLE"
    echo "Finished: $(date)"

    # Exit 0 so one failed sample does not mark the whole array task
    # as a SLURM failure. The failure is recorded in plasmer_status.tsv.
    exit 0
fi

# ------------------------------------------------------------
# Collect predicted plasmid FASTA
# ------------------------------------------------------------

if [ -f "$PRED_FASTA" ]; then
    cp "$PRED_FASTA" "$FINAL_FASTA"
    STATUS="OK"
else
    echo "WARNING: Predicted plasmid FASTA was not created:"
    echo "$PRED_FASTA"
    echo "Creating an empty final FASTA."

    touch "$FINAL_FASTA"
    STATUS="NO_PRED_FASTA"
fi

# Count FASTA records safely, including empty FASTA files
PLASMID_COUNT=$(
    awk '/^>/{count++} END{print count+0}' "$FINAL_FASTA"
)

# ------------------------------------------------------------
# Write completion status
# ------------------------------------------------------------

{
    echo -e "Sample\tPlatform\tSpecies\tStatus\tExit_Code\tFASTA\tPlasmid_Count\tFinal_FASTA"
    echo -e "$SAMPLE\t$PLATFORM\t$SPECIES\t$STATUS\t0\t$FASTA\t$PLASMID_COUNT\t$FINAL_FASTA"
} > "$STATUS_FILE"

echo ""
echo "Plasmids found: $PLASMID_COUNT"
echo "Status: $STATUS"
echo "Final FASTA: $FINAL_FASTA"
echo "Status file: $STATUS_FILE"
echo "DONE: $SAMPLE"
echo "Finished: $(date)"
