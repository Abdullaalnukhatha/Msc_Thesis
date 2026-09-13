#!/bin/bash
#SBATCH --job-name=rfplasmid_array
#SBATCH --time=04:00:00
#SBATCH --mem=24G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-120%20
#SBATCH --output=/rds/projects/e/elhamsak-mbru-amr/tempdata/logs/rfplasmid_%A_%a.out
#SBATCH --error=/rds/projects/e/elhamsak-mbru-amr/tempdata/logs/rfplasmid_%A_%a.err

set -euo pipefail

# ============================================================
# RFPlasmid array
#
# Expected manifest columns:
# SampleID    Platform    Species    Group    FASTA
#
# Expected input structure:
# results/Best_Assemblies/
# ├── ONT/
# │   ├── ecoli/
# │   └── kpneumoniae/
# └── PacBio/
#     ├── ecoli/
#     └── kpneumoniae/
#
# FASTA files are stored directly inside each species directory.
# ============================================================

PROJ="/rds/projects/e/elhamsak-mbru-amr"

BEST_ROOT="$PROJ/results/Best_Assemblies"
MANIFEST="$BEST_ROOT/master_manifest.tsv"

ENV_DIR="$PROJ/env/rfplasmid_env"
ENV_PYTHON="$ENV_DIR/bin/python"
ENV_RFPLASMID="$ENV_DIR/bin/rfplasmid"
RF_DB="$ENV_DIR/lib/python3.8/site-packages/RFPlasmid"

OUT_ROOT="$PROJ/results/plasmid_benchmark/RFPlasmid"
TEMP_ROOT="$PROJ/tempdata"

THREADS="${SLURM_CPUS_PER_TASK:-8}"
RFPLASMID_SPECIES="Enterobacteriaceae"

# ----------------------------
# Validate array context
# ----------------------------

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: SLURM_ARRAY_TASK_ID is not defined."
    echo "Submit this script using sbatch."
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

export TERM=xterm

# Separate temporary directory for each array task.
export TMPDIR="$TEMP_ROOT/tmp/rfplasmid_${ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export TMP="$TMPDIR"
export TEMP="$TMPDIR"

mkdir -p "$TMPDIR"

cleanup() {
    rm -rf "$TMPDIR"
}

trap cleanup EXIT

# ----------------------------
# Run information
# ----------------------------

echo "========================================="
echo "Running RFPlasmid array"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "SLURM job ID: ${SLURM_JOB_ID:-NA}"
echo "Array job ID: ${SLURM_ARRAY_JOB_ID:-NA}"
echo "Array task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Project: $PROJ"
echo "Best assemblies root: $BEST_ROOT"
echo "Manifest: $MANIFEST"
echo "Output root: $OUT_ROOT"
echo "RFPlasmid DB/package dir: $RF_DB"
echo "Temporary directory: $TMPDIR"
echo "Threads: $THREADS"
echo "RFPlasmid species model: $RFPLASMID_SPECIES"
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

export PATH="$ENV_DIR/bin:$PATH"

echo
echo "Conda environment:"
echo "CONDA_PREFIX=${CONDA_PREFIX:-NA}"
echo "which python: $(command -v python)"
echo "which rfplasmid: $(command -v rfplasmid)"
echo "ENV_PYTHON: $ENV_PYTHON"
echo "ENV_RFPLASMID: $ENV_RFPLASMID"

"$ENV_PYTHON" -V
"$ENV_PYTHON" -c "import pandas as pd; print('pandas', pd.__version__)"

"$ENV_RFPLASMID" --help | head || true

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

if [[ ! -x "$ENV_PYTHON" ]]; then
    echo "ERROR: Environment Python not executable:"
    echo "$ENV_PYTHON"
    exit 1
fi

if [[ ! -x "$ENV_RFPLASMID" ]]; then
    echo "ERROR: RFPlasmid executable not found:"
    echo "$ENV_RFPLASMID"
    exit 1
fi

if [[ ! -d "$RF_DB" ]]; then
    echo "ERROR: RFPlasmid package/database directory not found:"
    echo "$RF_DB"
    exit 1
fi

if [[ ! -f "$RF_DB/plasmiddb_cge.dmnd" ]]; then
    echo "ERROR: RFPlasmid database not found:"
    echo "$RF_DB/plasmiddb_cge.dmnd"
    exit 1
fi

if [[ ! -f "$RF_DB/plasmiddb_total.dmnd" ]]; then
    echo "ERROR: RFPlasmid database not found:"
    echo "$RF_DB/plasmiddb_total.dmnd"
    exit 1
fi

# ----------------------------
# Read manifest
# ----------------------------
# The first line is the header and must be skipped.
# Fields are read as tab-separated values.
# Windows carriage-return characters are removed.

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
    echo "ERROR: No records found in manifest."
    exit 1
fi

echo
echo "Total manifest records: $TOTAL"

if (( SLURM_ARRAY_TASK_ID > TOTAL )); then
    echo "Array task $SLURM_ARRAY_TASK_ID is greater than total records $TOTAL."
    echo "Exiting without error."
    exit 0
fi

RECORD="${RECORDS[$((SLURM_ARRAY_TASK_ID - 1))]}"

IFS=$'\t' read -r \
    SAMPLE \
    PLATFORM \
    SPECIES \
    GROUP \
    MANIFEST_FASTA \
    <<< "$RECORD"

# Remove any remaining carriage returns.
SAMPLE="${SAMPLE//$'\r'/}"
PLATFORM="${PLATFORM//$'\r'/}"
SPECIES="${SPECIES//$'\r'/}"
GROUP="${GROUP//$'\r'/}"
MANIFEST_FASTA="${MANIFEST_FASTA//$'\r'/}"

echo
echo "Manifest record:"
echo "Sample: $SAMPLE"
echo "Platform: $PLATFORM"
echo "Species: $SPECIES"
echo "Group: $GROUP"
echo "Manifest FASTA: $MANIFEST_FASTA"

# ----------------------------
# Normalize platform directory
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
# Normalize species directory
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
# Use the manifest path directly if it is valid.
if [[ -n "$MANIFEST_FASTA" && -f "$MANIFEST_FASTA" ]]; then
    ASSEMBLY="$MANIFEST_FASTA"
fi

# Option 2:
# Treat the manifest value as relative to the project.
if [[ -z "$ASSEMBLY" && -n "$MANIFEST_FASTA" && -f "$PROJ/$MANIFEST_FASTA" ]]; then
    ASSEMBLY="$PROJ/$MANIFEST_FASTA"
fi

# Option 3:
# Try the manifest basename inside the flat Best_Assemblies directory.
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
# Search the flat directory using common sample-based FASTA names.
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
    echo "Searched directory: $ASSEMBLY_DIR"
    echo
    echo "First files in the expected directory:"
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
# Prepare output paths
# ----------------------------

SAMPLE_DIR="$OUT_ROOT/$PLATFORM_DIR/$SPECIES_DIR/$SAMPLE"
FINAL_FASTA="$SAMPLE_DIR/${SAMPLE}_plasmids.fasta"
STATUS_TSV="$SAMPLE_DIR/rfplasmid_status.tsv"

mkdir -p "$SAMPLE_DIR"

echo
echo "Removing previous RFPlasmid output for this sample, if present..."

rm -rf "$SAMPLE_DIR"/run*
rm -f "$FINAL_FASTA" "$STATUS_TSV"

# RFPlasmid expects a directory containing FASTA files.
TEMP_INPUT="$TMPDIR/rfplasmid_${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}_${SAMPLE}"

rm -rf "$TEMP_INPUT"
mkdir -p "$TEMP_INPUT"

cp "$ASSEMBLY" "$TEMP_INPUT/${SAMPLE}.fasta"

# ----------------------------
# Run RFPlasmid
# ----------------------------

echo
echo "Running RFPlasmid..."
echo "Working directory for RFPlasmid:"
echo "$RF_DB"

OLDPWD_RUN="$PWD"
cd "$RF_DB"

set +e

"$ENV_RFPLASMID" \
    --species "$RFPLASMID_SPECIES" \
    --input "$TEMP_INPUT" \
    --out "$SAMPLE_DIR/run" \
    --threads "$THREADS" \
    --jelly

RF_EXIT=$?

set -e

cd "$OLDPWD_RUN"

rm -rf "$TEMP_INPUT"

if [[ "$RF_EXIT" -ne 0 ]]; then
    echo "WARNING: RFPlasmid exited with non-zero status: $RF_EXIT"
fi

# ----------------------------
# Locate prediction table
# ----------------------------

echo
echo "Searching for RFPlasmid prediction table..."

PREDICTION_FILE="$(
    find "$SAMPLE_DIR" \
        -type f \
        \( \
            -name "prediction.csv" \
            -o -name "prediction_full.csv" \
        \) \
        -print |
        sort |
        head -n 1 || true
)"

echo "Prediction file: ${PREDICTION_FILE:-NONE}"

# ----------------------------
# Extract predicted plasmid contigs
# ----------------------------

if [[ -n "${PREDICTION_FILE:-}" && -f "$PREDICTION_FILE" ]]; then

    "$ENV_PYTHON" - "$PREDICTION_FILE" "$ASSEMBLY" "$FINAL_FASTA" <<'EOF'
import sys
import csv
from pathlib import Path

pred_file = Path(sys.argv[1])
fasta_in = Path(sys.argv[2])
fasta_out = Path(sys.argv[3])


def detect_delimiter(path):
    text = path.read_text(errors="ignore")[:4096]

    if text.count("\t") >= text.count(","):
        return "\t"

    return ","


delimiter = detect_delimiter(pred_file)
plasmid_ids = set()

with open(pred_file, "r", newline="", errors="ignore") as f:
    reader = csv.DictReader(f, delimiter=delimiter)

    if reader.fieldnames:
        for row in reader:
            lower = {
                str(k).strip().lower(): str(v).strip()
                for k, v in row.items()
                if k is not None and v is not None
            }

            contig = (
                lower.get("contigid")
                or lower.get("contig id")
                or lower.get("contig")
                or lower.get("sequence")
                or lower.get("seqid")
                or lower.get("id")
                or lower.get("name")
            )

            if not contig:
                continue

            score_val = None

            for key in [
                "votes plasmid",
                "votes_plasmid",
                "plasmid",
                "plasmid score",
                "plasmid_score",
                "probability plasmid",
                "probability_plasmid",
                "prediction",
                "label",
                "class",
            ]:
                if key in lower:
                    score_val = lower[key]
                    break

            if score_val is None:
                continue

            is_plasmid = False

            try:
                is_plasmid = float(score_val) > 0.5
            except ValueError:
                if score_val.lower() in {
                    "plasmid",
                    "p",
                    "true",
                    "yes",
                    "1",
                }:
                    is_plasmid = True

            if is_plasmid:
                plasmid_ids.add(contig)

with open(fasta_in) as fin, open(fasta_out, "w") as fout:
    write = False

    for line in fin:
        if line.startswith(">"):
            contig_id = line[1:].split()[0].strip()
            write = contig_id in plasmid_ids

        if write:
            fout.write(line)
EOF

else
    echo "WARNING: RFPlasmid prediction file not found."
    echo "Creating empty final FASTA."
    : > "$FINAL_FASTA"
fi

# ----------------------------
# Count plasmid contigs
# ----------------------------

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

# ----------------------------
# Set status
# ----------------------------

if [[ "$RF_EXIT" -ne 0 ]]; then
    STATUS="WARNING"
    REASON="rfplasmid_nonzero_exit"

elif [[ "$PLASMID_COUNT" -gt 0 ]]; then
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
Sample	Platform	Species	Status	Reason	Plasmid_Count	Final_FASTA
$SAMPLE	$PLATFORM_DIR	$SPECIES_DIR	$STATUS	$REASON	$PLASMID_COUNT	$FINAL_FASTA
EOF

echo
echo "========================================="
echo "RFPlasmid finished"
echo "Sample: $SAMPLE"
echo "Platform: $PLATFORM_DIR"
echo "Species: $SPECIES_DIR"
echo "Input assembly: $ASSEMBLY"
echo "RFPlasmid exit code: $RF_EXIT"
echo "Plasmid count: $PLASMID_COUNT"
echo "Final FASTA: $FINAL_FASTA"
echo
echo "Status:"
cat "$STATUS_TSV"
echo
echo "Finished: $(date)"
echo "========================================="
