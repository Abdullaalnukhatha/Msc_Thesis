#!/bin/bash
#SBATCH --job-name=plassembler_array
#SBATCH --partition=emerald-shared
#SBATCH --time=12:00:00
#SBATCH --mem=48G
#SBATCH --cpus-per-task=8
#SBATCH --array=1-120%20
#SBATCH --output=/rds/projects/e/elhamsak-mbru-amr/tempdata/plassembler_%A_%a.out
#SBATCH --error=/rds/projects/e/elhamsak-mbru-amr/tempdata/plassembler_%A_%a.err

set -euo pipefail

export LC_ALL=C
export TERM="${TERM:-dumb}"

PROJ="/rds/projects/e/elhamsak-mbru-amr"

MANIFEST="$PROJ/results/Best_Assemblies/master_manifest.tsv"
READS_ROOT="$PROJ/results/ProcessedFASTQ"

# Selected assembly FASTA files and Flye assembly_info files
BEST_ASSEMBLIES_ROOT="$PROJ/results/Best_Assemblies"

DB_DIR="$PROJ/databases/PLSDB"
OUTPUT_ROOT="$PROJ/results/plasmid_benchmark/Plassembler"
ENV_DIR="$PROJ/env/plassembler_env"

THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$PROJ/tempdata"
mkdir -p "$OUTPUT_ROOT"

echo "========================================="
echo "Plassembler array job"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Started: $(date)"
echo "Host: $(hostname)"
echo "Manifest: $MANIFEST"
echo "Reads root: $READS_ROOT"
echo "Best assemblies root: $BEST_ASSEMBLIES_ROOT"
echo "Database: $DB_DIR"
echo "Output root: $OUTPUT_ROOT"
echo "Environment: $ENV_DIR"
echo "Threads: $THREADS"
echo "========================================="

module purge
module load bluebear
module load bear-apps/2024a/live
module load Miniforge3/25.3.0-3

source /rds/homes/a/axa2273/miniforge3/etc/profile.d/conda.sh

# Activate the environment through PATH and use conda run below
export CONDA_PREFIX="$ENV_DIR"
export PATH="$ENV_DIR/bin:/usr/bin:/bin:$PATH"
hash -r

# ============================================================
# Validate required paths
# ============================================================

if [ ! -f "$MANIFEST" ]; then
    echo "ERROR: Manifest not found: $MANIFEST"
    exit 1
fi

if [ ! -d "$READS_ROOT" ]; then
    echo "ERROR: Reads root not found: $READS_ROOT"
    exit 1
fi

if [ ! -d "$BEST_ASSEMBLIES_ROOT" ]; then
    echo "ERROR: Best assemblies root not found: $BEST_ASSEMBLIES_ROOT"
    exit 1
fi

if [ ! -d "$ENV_DIR" ]; then
    echo "ERROR: Environment not found: $ENV_DIR"
    exit 1
fi

if [ ! -d "$DB_DIR" ]; then
    echo "ERROR: Plassembler database not found: $DB_DIR"
    exit 1
fi

# ============================================================
# Dependency checks
# ============================================================

echo ""
echo "Dependency check using conda run:"

conda run -p "$ENV_DIR" which plassembler
conda run -p "$ENV_DIR" which unicycler
conda run -p "$ENV_DIR" which flye
conda run -p "$ENV_DIR" which raven
conda run -p "$ENV_DIR" which samtools
conda run -p "$ENV_DIR" which minimap2
conda run -p "$ENV_DIR" which mash
conda run -p "$ENV_DIR" which dnaapler

echo ""
echo "Version check:"

conda run -p "$ENV_DIR" plassembler --version || true
conda run -p "$ENV_DIR" unicycler --version || true
conda run -p "$ENV_DIR" flye --version || true
conda run -p "$ENV_DIR" raven --version || true

# ============================================================
# Read manifest row
# ============================================================

echo ""
echo "Reading manifest row..."

ROW=$(
    tail -n +2 "$MANIFEST" |
    sed -n "${SLURM_ARRAY_TASK_ID}p" |
    tr -d '\r'
)

if [ -z "$ROW" ]; then
    echo "ERROR: No row found for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

SAMPLE=$(printf '%s\n' "$ROW" | cut -f1)
PLATFORM=$(printf '%s\n' "$ROW" | cut -f2)
SPECIES=$(printf '%s\n' "$ROW" | cut -f3)
GROUP=$(printf '%s\n' "$ROW" | cut -f4)
ASSEMBLY_FASTA=$(printf '%s\n' "$ROW" | cut -f5)

READS="$READS_ROOT/$PLATFORM/$SPECIES/${SAMPLE}_processed.fastq.gz"

# Expected example:
# /results/Best_Assemblies/PacBio/ecoli/SRR7647824_assembly_info.txt
FLYE_INFO="$BEST_ASSEMBLIES_ROOT/$PLATFORM/$SPECIES/${SAMPLE}_assembly_info.txt"

SAMPLE_OUT="$OUTPUT_ROOT/$PLATFORM/$SPECIES/$SAMPLE"
RUN_OUT="$SAMPLE_OUT/run"
FINAL_FASTA="$SAMPLE_OUT/${SAMPLE}_plasmids.fasta"
STATUS_FILE="$SAMPLE_OUT/plassembler_status.tsv"

mkdir -p "$SAMPLE_OUT"

echo ""
echo "Parsed values:"
echo "  SAMPLE        = [$SAMPLE]"
echo "  PLATFORM      = [$PLATFORM]"
echo "  SPECIES       = [$SPECIES]"
echo "  GROUP         = [$GROUP]"
echo "  READS         = [$READS]"
echo "  ASSEMBLY      = [$ASSEMBLY_FASTA]"
echo "  FLYE_INFO     = [$FLYE_INFO]"
echo "  OUTDIR        = [$SAMPLE_OUT]"

# ============================================================
# Skip samples that completed successfully
# ============================================================
#
# The status file is used instead of checking only whether the
# final FASTA is non-empty.
#
# A valid completed sample with zero plasmids has:
#   Status = OK
#   Reason = no_plasmids_detected
#   An empty final FASTA
#
# Failed or incomplete samples are rerun.
# ============================================================

if [ -s "$STATUS_FILE" ]; then
    EXISTING_STATUS=$(
        awk -F'\t' 'NR == 2 {gsub(/\r/, "", $4); print $4}' \
            "$STATUS_FILE"
    )

    EXISTING_REASON=$(
        awk -F'\t' 'NR == 2 {gsub(/\r/, "", $5); print $5}' \
            "$STATUS_FILE"
    )

    EXISTING_COUNT=$(
        awk -F'\t' 'NR == 2 {gsub(/\r/, "", $6); print $6}' \
            "$STATUS_FILE"
    )

    if [ "$EXISTING_STATUS" = "OK" ]; then
        echo ""
        echo "SKIPPING: Sample already completed successfully."
        echo "  Sample:        $SAMPLE"
        echo "  Status:        $EXISTING_STATUS"
        echo "  Reason:        ${EXISTING_REASON:-unknown}"
        echo "  Plasmid count: ${EXISTING_COUNT:-unknown}"
        echo "  Final FASTA:   $FINAL_FASTA"
        echo "  Status file:   $STATUS_FILE"
        echo "Finished: $(date)"
        exit 0
    fi

    echo ""
    echo "Previous result was not successful and will be reprocessed."
    echo "  Previous status: ${EXISTING_STATUS:-unknown}"
    echo "  Previous reason: ${EXISTING_REASON:-unknown}"
fi

# ============================================================
# Validate input reads
# ============================================================

if [ ! -s "$READS" ]; then
    echo "ERROR: Reads file missing or empty: $READS"

    printf \
        "Sample\tPlatform\tSpecies\tStatus\tReason\tPlasmid_Count\tFinal_FASTA\n" \
        > "$STATUS_FILE"

    printf \
        "%s\t%s\t%s\tFAILED\tmissing_reads\t0\t%s\n" \
        "$SAMPLE" \
        "$PLATFORM" \
        "$SPECIES" \
        "$FINAL_FASTA" \
        >> "$STATUS_FILE"

    # Create or empty the final FASTA
    : > "$FINAL_FASTA"

    echo "DONE_WITH_FAILURE_RECORDED: $SAMPLE"
    exit 0
fi

# ============================================================
# Validate assembly FASTA
# ============================================================

if [ ! -s "$ASSEMBLY_FASTA" ]; then
    echo "ERROR: Assembly FASTA missing or empty: $ASSEMBLY_FASTA"

    printf \
        "Sample\tPlatform\tSpecies\tStatus\tReason\tPlasmid_Count\tFinal_FASTA\n" \
        > "$STATUS_FILE"

    printf \
        "%s\t%s\t%s\tFAILED\tmissing_assembly\t0\t%s\n" \
        "$SAMPLE" \
        "$PLATFORM" \
        "$SPECIES" \
        "$FINAL_FASTA" \
        >> "$STATUS_FILE"

    : > "$FINAL_FASTA"

    echo "DONE_WITH_FAILURE_RECORDED: $SAMPLE"
    exit 0
fi

# ============================================================
# Validate Flye assembly_info.txt
# ============================================================

if [ ! -s "$FLYE_INFO" ]; then
    echo "ERROR: Flye assembly_info.txt missing or empty: $FLYE_INFO"

    printf \
        "Sample\tPlatform\tSpecies\tStatus\tReason\tPlasmid_Count\tFinal_FASTA\n" \
        > "$STATUS_FILE"

    printf \
        "%s\t%s\t%s\tFAILED\tmissing_flye_info\t0\t%s\n" \
        "$SAMPLE" \
        "$PLATFORM" \
        "$SPECIES" \
        "$FINAL_FASTA" \
        >> "$STATUS_FILE"

    : > "$FINAL_FASTA"

    echo "DONE_WITH_FAILURE_RECORDED: $SAMPLE"
    exit 0
fi

# ============================================================
# Validate compressed FASTQ
# ============================================================

if ! gzip -t "$READS" 2>/dev/null; then
    echo "ERROR: Reads file is corrupted: $READS"

    printf \
        "Sample\tPlatform\tSpecies\tStatus\tReason\tPlasmid_Count\tFinal_FASTA\n" \
        > "$STATUS_FILE"

    printf \
        "%s\t%s\t%s\tFAILED\tcorrupted_reads\t0\t%s\n" \
        "$SAMPLE" \
        "$PLATFORM" \
        "$SPECIES" \
        "$FINAL_FASTA" \
        >> "$STATUS_FILE"

    : > "$FINAL_FASTA"

    echo "DONE_WITH_FAILURE_RECORDED: $SAMPLE"
    exit 0
fi

# ============================================================
# Choose platform-specific Plassembler options
# ============================================================

case "$PLATFORM" in
    ONT)
        PLATFORM_ARGS=(--raw_flag)
        ;;

    PacBio)
        PLATFORM_ARGS=(--pacbio_model pacbio-hifi)
        ;;

    *)
        echo "ERROR: Unknown platform: $PLATFORM"

        printf \
            "Sample\tPlatform\tSpecies\tStatus\tReason\tPlasmid_Count\tFinal_FASTA\n" \
            > "$STATUS_FILE"

        printf \
            "%s\t%s\t%s\tFAILED\tunknown_platform\t0\t%s\n" \
            "$SAMPLE" \
            "$PLATFORM" \
            "$SPECIES" \
            "$FINAL_FASTA" \
            >> "$STATUS_FILE"

        : > "$FINAL_FASTA"

        echo "DONE_WITH_FAILURE_RECORDED: $SAMPLE"
        exit 0
        ;;
esac

# ============================================================
# Run Plassembler
# ============================================================

echo ""
echo "Running Plassembler..."
echo "Command platform arguments: ${PLATFORM_ARGS[*]}"

# Remove only the incomplete or failed previous run
rm -rf "$RUN_OUT"
rm -f "$FINAL_FASTA"

set +e

conda run -p "$ENV_DIR" plassembler long \
    -d "$DB_DIR" \
    -l "$READS" \
    -o "$RUN_OUT" \
    --skip_qc \
    -t "$THREADS" \
    "${PLATFORM_ARGS[@]}" \
    --flye_assembly "$ASSEMBLY_FASTA" \
    --flye_info "$FLYE_INFO" \
    --force

PLASSEMBLER_EXIT=$?

set -e

# ============================================================
# Check Plassembler output
# ============================================================

echo ""
echo "Checking Plassembler output..."

CANDIDATE_1="$RUN_OUT/plassembler_plasmids.fasta"
CANDIDATE_2="$RUN_OUT/plasmids.fasta"
CANDIDATE_3="$RUN_OUT/${SAMPLE}_plasmids.fasta"

if [ "$PLASSEMBLER_EXIT" -ne 0 ]; then
    echo "WARNING: Plassembler failed for $SAMPLE"
    echo "Exit code: $PLASSEMBLER_EXIT"

    : > "$FINAL_FASTA"

    printf \
        "Sample\tPlatform\tSpecies\tStatus\tReason\tPlasmid_Count\tFinal_FASTA\n" \
        > "$STATUS_FILE"

    printf \
        "%s\t%s\t%s\tFAILED\tplassembler_exit_%s\t0\t%s\n" \
        "$SAMPLE" \
        "$PLATFORM" \
        "$SPECIES" \
        "$PLASSEMBLER_EXIT" \
        "$FINAL_FASTA" \
        >> "$STATUS_FILE"

    echo "Plasmids found: 0"
    echo "Final FASTA: $FINAL_FASTA"
    echo "Status file: $STATUS_FILE"
    echo "DONE_WITH_FAILURE_RECORDED: $SAMPLE"
    echo "Finished: $(date)"

    exit 0
fi

# ============================================================
# Copy the detected plasmid FASTA
# ============================================================

if [ -s "$CANDIDATE_1" ]; then
    echo "Using output: $CANDIDATE_1"
    cp "$CANDIDATE_1" "$FINAL_FASTA"

elif [ -s "$CANDIDATE_2" ]; then
    echo "Using output: $CANDIDATE_2"
    cp "$CANDIDATE_2" "$FINAL_FASTA"

elif [ -s "$CANDIDATE_3" ]; then
    echo "Using output: $CANDIDATE_3"
    cp "$CANDIDATE_3" "$FINAL_FASTA"

else
    echo "No plasmid FASTA output found."
    echo "Recording the completed sample as zero plasmids."

    : > "$FINAL_FASTA"
fi

# Count FASTA records safely
PLASMID_COUNT=$(
    awk '/^>/ {count++} END {print count + 0}' "$FINAL_FASTA"
)

if [ "$PLASMID_COUNT" -gt 0 ]; then
    STATUS="OK"
    REASON="plasmids_detected"
else
    STATUS="OK"
    REASON="no_plasmids_detected"
fi

# ============================================================
# Write completion status
# ============================================================

printf \
    "Sample\tPlatform\tSpecies\tStatus\tReason\tPlasmid_Count\tFinal_FASTA\n" \
    > "$STATUS_FILE"

printf \
    "%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "$SAMPLE" \
    "$PLATFORM" \
    "$SPECIES" \
    "$STATUS" \
    "$REASON" \
    "$PLASMID_COUNT" \
    "$FINAL_FASTA" \
    >> "$STATUS_FILE"

echo ""
echo "Plasmids found: $PLASMID_COUNT"
echo "Final FASTA: $FINAL_FASTA"
echo "Status file: $STATUS_FILE"
echo "DONE: $SAMPLE"
echo "Finished: $(date)"
