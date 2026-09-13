#!/bin/bash
#SBATCH --job-name=platon_array
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-120%20
#SBATCH --output=/rds/projects/e/elhamsak-mbru-amr/tempdata/logs/platon_%A_%a.out
#SBATCH --error=/rds/projects/e/elhamsak-mbru-amr/tempdata/logs/platon_%A_%a.err

set -euo pipefail

# ============================================================
# Platon array
#
# Manifest columns:
# SampleID    Platform    Species    Group    FASTA
#
# Input structure:
# results/Best_Assemblies/
# ├── ONT/
# │   ├── ecoli/
# │   └── kpneumoniae/
# └── PacBio/
#     ├── ecoli/
#     └── kpneumoniae/
#
# FASTA files are stored directly inside each species directory.
# There is no separate input folder per sample.
# ============================================================

PROJ="/rds/projects/e/elhamsak-mbru-amr"

BEST_ROOT="$PROJ/results/Best_Assemblies"
MANIFEST="$BEST_ROOT/master_manifest.tsv"

ENV_DIR="$PROJ/env/plasmidEC_platon"
OUT_ROOT="$PROJ/results/plasmid_benchmark/Platon"
TEMP_ROOT="$PROJ/tempdata"

DB="$PROJ/databases/PlasmerDB/platon_db"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

# ----------------------------
# Validate SLURM array context
# ----------------------------

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not defined."
    echo "Submit this script with sbatch."
    exit 1
fi

ARRAY_JOB_ID="${SLURM_ARRAY_JOB_ID:-${SLURM_JOB_ID:-manual}}"

# ----------------------------
# Create directories
# ----------------------------

mkdir -p \
    "$OUT_ROOT" \
    "$TEMP_ROOT/logs" \
    "$TEMP_ROOT/tmp"

# Give every array task its own temporary directory.
export TERM=xterm
export TMPDIR="$TEMP_ROOT/tmp/platon_${ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export TEMP="$TMPDIR"
export TMP="$TMPDIR"

mkdir -p "$TMPDIR"

cleanup() {
    rm -rf "$TMPDIR"
}

trap cleanup EXIT

# ----------------------------
# Run information
# ----------------------------

echo "========================================="
echo "Running Platon array"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "SLURM job ID: ${SLURM_JOB_ID:-NA}"
echo "Array job ID: ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Project: $PROJ"
echo "Best assemblies root: $BEST_ROOT"
echo "Manifest: $MANIFEST"
echo "Output root: $OUT_ROOT"
echo "Temporary directory: $TMPDIR"
echo "Threads: $THREADS"
echo "Database: $DB"
echo "========================================="

# ----------------------------
# Load environment
# ----------------------------

module purge
module load bluebear
module load bear-apps/2024a/live
module load Miniforge3/25.3.0-3

source /rds/homes/a/axa2273/miniforge3/etc/profile.d/conda.sh
conda activate "$ENV_DIR"

echo
echo "Conda environment:"
echo "Environment: ${CONDA_PREFIX:-NA}"
echo "Platon executable: $(command -v platon)"
platon --version || true
echo

# ----------------------------
# Validate required inputs
# ----------------------------

if [[ ! -d "$BEST_ROOT" ]]; then
    echo "ERROR: Best assemblies directory not found:"
    echo "$BEST_ROOT"
    exit 1
fi

if [[ ! -f "$MANIFEST" ]]; then
    echo "ERROR: Manifest not found:"
    echo "$MANIFEST"
    exit 1
fi

if [[ ! -d "$DB" ]]; then
    echo "ERROR: Platon database not found:"
    echo "$DB"
    exit 1
fi

if ! command -v platon >/dev/null 2>&1; then
    echo "ERROR: Platon was not found in the active environment."
    exit 1
fi

# ----------------------------
# Read manifest
# ----------------------------
# Skip line 1 because it is the header.
# Require at least five tab-separated columns.
# Remove Windows carriage-return characters.

mapfile -t RECORDS < <(
    awk -F'\t' '
        BEGIN {
            OFS = "\t"
        }

        NR == 1 {
            next
        }

        NF >= 5 && $1 != "" {
            gsub(/\r/, "", $0)
            print $1, $2, $3, $4, $5
        }
    ' "$MANIFEST"
)

TOTAL="${#RECORDS[@]}"

if [[ "$TOTAL" -eq 0 ]]; then
    echo "ERROR: No sample records were found in the manifest."
    exit 1
fi

echo "Total manifest records: $TOTAL"

if (( SLURM_ARRAY_TASK_ID > TOTAL )); then
    echo "Array task $SLURM_ARRAY_TASK_ID exceeds total records $TOTAL."
    echo "Exiting without error."
    exit 0
fi

# Bash arrays start at index 0, while SLURM tasks start at 1.
RECORD="${RECORDS[$((SLURM_ARRAY_TASK_ID - 1))]}"

IFS=$'\t' read -r \
    SAMPLE \
    PLATFORM \
    SPECIES \
    GROUP \
    MANIFEST_FASTA \
    <<< "$RECORD"

# Remove any remaining carriage-return characters.
SAMPLE="${SAMPLE//$'\r'/}"
PLATFORM="${PLATFORM//$'\r'/}"
SPECIES="${SPECIES//$'\r'/}"
GROUP="${GROUP//$'\r'/}"
MANIFEST_FASTA="${MANIFEST_FASTA//$'\r'/}"

if [[ -z "$SAMPLE" || -z "$PLATFORM" || -z "$SPECIES" ]]; then
    echo "ERROR: Manifest record contains an empty required field."
    echo "Record: $RECORD"
    exit 1
fi

echo
echo "Manifest record:"
echo "Sample: $SAMPLE"
echo "Platform: $PLATFORM"
echo "Species: $SPECIES"
echo "Group: $GROUP"
echo "Manifest FASTA: $MANIFEST_FASTA"

# ----------------------------
# Normalize platform
# ----------------------------

case "${PLATFORM,,}" in
    ont|nanopore|oxford_nanopore|oxfordnanopore)
        PLATFORM_DIR="ONT"
        ;;

    pacbio|pacbio_hifi|pacbiohifi|hifi)
        PLATFORM_DIR="PacBio"
        ;;

    *)
        PLATFORM_DIR="$PLATFORM"
        ;;
esac

# ----------------------------
# Normalize species
# ----------------------------

case "${SPECIES,,}" in
    ecoli|e_coli|escherichia_coli|escherichia-coli)
        SPECIES_DIR="ecoli"
        ;;

    kpneumoniae|k_pneumoniae|pneumoniae|klebsiella_pneumoniae|klebsiella-pneumoniae)
        SPECIES_DIR="kpneumoniae"
        ;;

    *)
        SPECIES_DIR="$SPECIES"
        ;;
esac

ASSEMBLY_DIR="$BEST_ROOT/$PLATFORM_DIR/$SPECIES_DIR"

echo
echo "Resolved platform directory: $PLATFORM_DIR"
echo "Resolved species directory: $SPECIES_DIR"
echo "Flat assembly directory: $ASSEMBLY_DIR"

if [[ ! -d "$ASSEMBLY_DIR" ]]; then
    echo "ERROR: Expected assembly directory does not exist:"
    echo "$ASSEMBLY_DIR"
    exit 1
fi

# ----------------------------
# Resolve assembly FASTA
# ----------------------------

ASSEMBLY=""

# Option 1:
# Use the manifest FASTA directly if it is a valid absolute path.
if [[ -n "$MANIFEST_FASTA" && -f "$MANIFEST_FASTA" ]]; then
    ASSEMBLY="$MANIFEST_FASTA"
fi

# Option 2:
# Try the manifest value as a path relative to the project.
if [[ -z "$ASSEMBLY" && -n "$MANIFEST_FASTA" && -f "$PROJ/$MANIFEST_FASTA" ]]; then
    ASSEMBLY="$PROJ/$MANIFEST_FASTA"
fi

# Option 3:
# Try the manifest FASTA basename directly in the flat assembly directory.
if [[ -z "$ASSEMBLY" && -n "$MANIFEST_FASTA" ]]; then
    FASTA_BASENAME="$(basename "$MANIFEST_FASTA")"

    if [[ -f "$ASSEMBLY_DIR/$FASTA_BASENAME" ]]; then
        ASSEMBLY="$ASSEMBLY_DIR/$FASTA_BASENAME"
    fi
fi

# Option 4:
# Search the flat directory for the exact manifest basename.
if [[ -z "$ASSEMBLY" && -n "$MANIFEST_FASTA" ]]; then
    FASTA_BASENAME="$(basename "$MANIFEST_FASTA")"

    mapfile -t EXACT_MATCHES < <(
        find "$ASSEMBLY_DIR" \
            -maxdepth 1 \
            -type f \
            -name "$FASTA_BASENAME" \
            -print
    )

    if [[ "${#EXACT_MATCHES[@]}" -eq 1 ]]; then
        ASSEMBLY="${EXACT_MATCHES[0]}"

    elif [[ "${#EXACT_MATCHES[@]}" -gt 1 ]]; then
        echo "ERROR: Multiple FASTA files matched the manifest basename:"
        printf '  %s\n' "${EXACT_MATCHES[@]}"
        exit 1
    fi
fi

# Option 5:
# Search the flat directory using common SampleID-based filenames.
if [[ -z "$ASSEMBLY" ]]; then
    mapfile -t SAMPLE_MATCHES < <(
        find "$ASSEMBLY_DIR" \
            -maxdepth 1 \
            -type f \
            \( \
                -iname "${SAMPLE}.fasta" \
                -o -iname "${SAMPLE}.fa" \
                -o -iname "${SAMPLE}.fna" \
                -o -iname "${SAMPLE}_assembly.fasta" \
                -o -iname "${SAMPLE}_assembly.fa" \
                -o -iname "${SAMPLE}_assembly.fna" \
                -o -iname "${SAMPLE}_contigs.fasta" \
                -o -iname "${SAMPLE}_contigs.fa" \
                -o -iname "${SAMPLE}_contigs.fna" \
                -o -iname "${SAMPLE}*.fasta" \
                -o -iname "${SAMPLE}*.fa" \
                -o -iname "${SAMPLE}*.fna" \
            \) \
            -print
    )

    if [[ "${#SAMPLE_MATCHES[@]}" -eq 1 ]]; then
        ASSEMBLY="${SAMPLE_MATCHES[0]}"

    elif [[ "${#SAMPLE_MATCHES[@]}" -gt 1 ]]; then
        echo "ERROR: Multiple assembly FASTA files matched sample $SAMPLE:"
        printf '  %s\n' "${SAMPLE_MATCHES[@]}"
        echo
        echo "Update the manifest FASTA column with the exact required path."
        exit 1
    fi
fi

if [[ -z "$ASSEMBLY" || ! -f "$ASSEMBLY" ]]; then
    echo "ERROR: Could not locate assembly FASTA."
    echo
    echo "Sample: $SAMPLE"
    echo "Manifest FASTA: $MANIFEST_FASTA"
    echo "Searched flat directory: $ASSEMBLY_DIR"
    echo
    echo "First files present in the directory:"
    find "$ASSEMBLY_DIR" \
        -maxdepth 1 \
        -type f \
        -printf '%f\n' |
        head -30
    exit 1
fi

if [[ ! -s "$ASSEMBLY" ]]; then
    echo "ERROR: Assembly FASTA is empty:"
    echo "$ASSEMBLY"
    exit 1
fi

echo "Resolved assembly FASTA: $ASSEMBLY"

# ----------------------------
# Prepare output directory
# ----------------------------

SAMPLE_DIR="$OUT_ROOT/$PLATFORM_DIR/$SPECIES_DIR/$SAMPLE"
RUN_DIR="$SAMPLE_DIR/run"

FINAL_FASTA="$SAMPLE_DIR/${SAMPLE}_plasmids.fasta"
STATUS_TSV="$SAMPLE_DIR/platon_status.tsv"

mkdir -p "$SAMPLE_DIR"

echo
echo "Removing previous Platon output for this sample, if present..."

rm -rf "$RUN_DIR"
rm -f "$FINAL_FASTA" "$STATUS_TSV"

mkdir -p "$RUN_DIR"

# ----------------------------
# Run Platon
# ----------------------------

echo
echo "Running Platon..."

platon \
    --db "$DB" \
    --output "$RUN_DIR" \
    --prefix "$SAMPLE" \
    --mode accuracy \
    --threads "$THREADS" \
    "$ASSEMBLY"

# ----------------------------
# Locate Platon output
# ----------------------------

PLATON_FASTA="$RUN_DIR/${SAMPLE}.plasmid.fasta"

if [[ -s "$PLATON_FASTA" ]]; then
    cp "$PLATON_FASTA" "$FINAL_FASTA"
else
    echo "No non-empty Platon plasmid FASTA was produced."
    : > "$FINAL_FASTA"
fi

PLASMID_COUNT="$(
    awk '
        /^>/ {
            count++
        }

        END {
            print count + 0
        }
    ' "$FINAL_FASTA"
)"

if (( PLASMID_COUNT > 0 )); then
    STATUS="OK"
    REASON="plasmids_detected"
else
    STATUS="OK"
    REASON="no_plasmids_detected"
fi

# ----------------------------
# Write status file
# ----------------------------

cat > "$STATUS_TSV" <<EOF
Sample	Platform	Species	Group	Status	Reason	Plasmid_Count	Input_FASTA	Final_FASTA
$SAMPLE	$PLATFORM_DIR	$SPECIES_DIR	$GROUP	$STATUS	$REASON	$PLASMID_COUNT	$ASSEMBLY	$FINAL_FASTA
EOF

echo
echo "========================================="
echo "Platon finished"
echo "Sample: $SAMPLE"
echo "Platform: $PLATFORM_DIR"
echo "Species: $SPECIES_DIR"
echo "Input FASTA: $ASSEMBLY"
echo "Plasmid count: $PLASMID_COUNT"
echo "Final FASTA: $FINAL_FASTA"
echo
echo "Status:"
cat "$STATUS_TSV"
echo
echo "Finished: $(date)"
echo "========================================="
