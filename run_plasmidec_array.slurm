#!/bin/bash
#SBATCH --job-name=plasmidec_array
#SBATCH --time=04:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-120%20
#SBATCH --output=/rds/projects/e/elhamsak-mbru-amr/tempdata/logs/plasmidec_%A_%a.out
#SBATCH --error=/rds/projects/e/elhamsak-mbru-amr/tempdata/logs/plasmidec_%A_%a.err

set -euo pipefail

PROJ="/rds/projects/e/elhamsak-mbru-amr"

MANIFEST="$PROJ/results/Best_Assemblies/master_manifest.tsv"
BEST_ASSEMBLIES_ROOT="$PROJ/results/Best_Assemblies"

PLASMIDEC_DIR="$PROJ/tools/plasmidEC"

CENTRIFUGE_ENV="$PROJ/env/plasmidEC_centrifuge"
PLATON_ENV="$PROJ/env/plasmidEC_platon"
R_ENV="$PROJ/env/plasmidEC_R"
RFPLASMID_ENV="$PROJ/env/rfplasmid_env"

OUT_ROOT="$PROJ/results/plasmid_benchmark/PlasmidEC"
TEMP_ROOT="$PROJ/tempdata"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$OUT_ROOT"
mkdir -p "$TEMP_ROOT/logs"
mkdir -p "$TEMP_ROOT/tmp"

export TERM="${TERM:-xterm}"
export TMPDIR="$TEMP_ROOT/tmp"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"

echo "========================================="
echo "Running PlasmidEC array"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "Array job ID: ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID:-NA}"
echo "Project: $PROJ"
echo "Manifest: $MANIFEST"
echo "Best assemblies root: $BEST_ASSEMBLIES_ROOT"
echo "Output root: $OUT_ROOT"
echo "PlasmidEC dir: $PLASMIDEC_DIR"
echo "Threads: $THREADS"
echo "========================================="

module purge
module load bluebear
module load bear-apps/2024a/live
module load Miniforge3/25.3.0-3

source /rds/homes/a/axa2273/miniforge3/etc/profile.d/conda.sh

if [ ! -f "$MANIFEST" ]; then
    echo "ERROR: Manifest not found:"
    echo "$MANIFEST"
    exit 1
fi

if [ ! -d "$BEST_ASSEMBLIES_ROOT" ]; then
    echo "ERROR: Best assemblies directory not found:"
    echo "$BEST_ASSEMBLIES_ROOT"
    exit 1
fi

if [ ! -f "$PLASMIDEC_DIR/plasmidEC.sh" ]; then
    echo "ERROR: plasmidEC.sh not found:"
    echo "$PLASMIDEC_DIR/plasmidEC.sh"
    exit 1
fi

echo ""
echo "Checking PlasmidEC environments..."

for ENV_PATH in \
    "$CENTRIFUGE_ENV" \
    "$PLATON_ENV" \
    "$R_ENV" \
    "$RFPLASMID_ENV"
do
    if [ ! -d "$ENV_PATH" ]; then
        echo "ERROR: Required PlasmidEC environment not found:"
        echo "$ENV_PATH"
        exit 1
    fi

    echo "Found: $ENV_PATH"
done

echo ""
echo "Checking PlasmidEC script:"
ls -lh "$PLASMIDEC_DIR/plasmidEC.sh"

echo ""
echo "Reading manifest..."

mapfile -t RECORDS < <(
    awk -F'\t' '
        BEGIN {
            OFS = "\t"
        }

        NR > 1 && NF >= 5 {
            gsub(/\r$/, "", $5)
            print $1, $2, $3, $4, $5
        }
    ' "$MANIFEST"
)

TOTAL="${#RECORDS[@]}"

if [ "$TOTAL" -eq 0 ]; then
    echo "ERROR: No sample records found in manifest."
    exit 1
fi

echo "Total manifest records: $TOTAL"

if [ "$SLURM_ARRAY_TASK_ID" -gt "$TOTAL" ]; then
    echo "Array task $SLURM_ARRAY_TASK_ID is greater than total records $TOTAL."
    echo "Exiting without error."
    exit 0
fi

RECORD="${RECORDS[$((SLURM_ARRAY_TASK_ID - 1))]}"

IFS=$'\t' read -r SAMPLE PLATFORM SPECIES GROUP ASSEMBLY <<< "$RECORD"

echo ""
echo "Manifest values:"
echo "Sample: $SAMPLE"
echo "Platform: $PLATFORM"
echo "Species: $SPECIES"
echo "Group: $GROUP"
echo "Assembly: $ASSEMBLY"

if [ -z "$SAMPLE" ] || \
   [ -z "$PLATFORM" ] || \
   [ -z "$SPECIES" ] || \
   [ -z "$GROUP" ] || \
   [ -z "$ASSEMBLY" ]; then

    echo "ERROR: One or more manifest fields are empty."
    exit 1
fi

case "$ASSEMBLY" in
    "$BEST_ASSEMBLIES_ROOT"/*)
        ;;
    *)
        echo "ERROR: Assembly path is outside Best_Assemblies:"
        echo "$ASSEMBLY"
        exit 1
        ;;
esac

if [ ! -s "$ASSEMBLY" ]; then
    echo "ERROR: Assembly FASTA is missing or empty:"
    echo "$ASSEMBLY"
    exit 1
fi

SAMPLE_DIR="$OUT_ROOT/$PLATFORM/$SPECIES/$SAMPLE"
RUN_DIR="$SAMPLE_DIR/run"
FINAL_FASTA="$SAMPLE_DIR/${SAMPLE}_plasmids.fasta"
STATUS_TSV="$SAMPLE_DIR/plasmidec_status.tsv"

mkdir -p "$SAMPLE_DIR"

echo ""
echo "Removing previous output for this sample, if present..."

rm -rf "$RUN_DIR"
rm -f "$FINAL_FASTA"
rm -f "$STATUS_TSV"

# Do not create RUN_DIR here.
# PlasmidEC creates its own output directory.

echo ""
echo "Running PlasmidEC..."

set +e

bash "$PLASMIDEC_DIR/plasmidEC.sh" \
    -i "$ASSEMBLY" \
    -o "$RUN_DIR" \
    -n "$SAMPLE" \
    -l 500 \
    -t "$THREADS" \
    -f

PEC_EXIT=$?

set -e

echo ""
echo "PlasmidEC exit code: $PEC_EXIT"

if [ "$PEC_EXIT" -ne 0 ]; then
    echo "WARNING: PlasmidEC exited with non-zero status."
fi

echo ""
echo "Searching for PlasmidEC plasmid FASTA..."

PLASMIDEC_FASTA=$(find "$RUN_DIR" \
    -type f \
    \( \
        -name "plasmid_predictions.fasta" -o \
        -name "*plasmid*.fasta" -o \
        -name "*plasmid*.fa" -o \
        -name "*plasmid*.fna" \
    \) \
    -size +0c \
    | sort \
    | head -n 1 \
    || true)

echo "PlasmidEC FASTA found: ${PLASMIDEC_FASTA:-NONE}"

if [ -n "${PLASMIDEC_FASTA:-}" ] && [ -s "$PLASMIDEC_FASTA" ]; then
    cp "$PLASMIDEC_FASTA" "$FINAL_FASTA"
else
    touch "$FINAL_FASTA"
fi

PLASMID_COUNT=$(awk '/^>/ {count++} END {print count+0}' "$FINAL_FASTA")

if [ "$PEC_EXIT" -ne 0 ]; then
    STATUS="WARNING"
    REASON="plasmidec_nonzero_exit_${PEC_EXIT}"

elif [ "$PLASMID_COUNT" -gt 0 ]; then
    STATUS="OK"
    REASON="plasmids_detected"

else
    STATUS="OK"
    REASON="no_plasmids_detected"
fi

cat > "$STATUS_TSV" <<EOF
Sample	Platform	Species	Status	Reason	Plasmid_Count	Final_FASTA
$SAMPLE	$PLATFORM	$SPECIES	$STATUS	$REASON	$PLASMID_COUNT	$FINAL_FASTA
EOF

echo ""
echo "========================================="
echo "PlasmidEC finished"
echo "Sample: $SAMPLE"
echo "Platform: $PLATFORM"
echo "Species: $SPECIES"
echo "Exit code: $PEC_EXIT"
echo "Plasmid count: $PLASMID_COUNT"
echo "Final FASTA: $FINAL_FASTA"
echo "Status file: $STATUS_TSV"
echo ""
cat "$STATUS_TSV"
echo ""
echo "Finished: $(date)"
echo "========================================="
