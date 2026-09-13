#!/bin/bash
# Align selected best assembly contigs against their full reference genomes
# Project: /rds/projects/e/elhamsak-mbru-amr

#SBATCH --job-name=align_assemblies
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=/rds/projects/e/elhamsak-mbru-amr/tempdata/align_assemblies_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-mbru-amr/tempdata/align_assemblies_%j.err

set -Eeuo pipefail
trap 'echo "ERROR at line $LINENO: $BASH_COMMAND" >&2' ERR

# =====================================================================
# Project directories
# =====================================================================
PROJECT_ROOT="/rds/projects/e/elhamsak-mbru-amr"
ASSEMBLY_DIR="$PROJECT_ROOT/results/Best_Assemblies"
REFERENCE_DIR="$PROJECT_ROOT/data/ReferenceFASTA"
OUTPUT_DIR="$PROJECT_ROOT/results/alignments"
MAPPING_FILE="$REFERENCE_DIR/reference_mapping.tsv"
TEMP_DIR="$PROJECT_ROOT/tempdata"
CONDA_ENV="$PROJECT_ROOT/envs/ref_valid"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$OUTPUT_DIR" "$TEMP_DIR"

# =====================================================================
# BlueBEAR environment
# =====================================================================
module purge
module load bluebear
module load bear-apps/2024a/live
module load Miniforge3/25.3.0-3

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

command -v minimap2 >/dev/null 2>&1 || {
    echo "ERROR: minimap2 was not found after activating $CONDA_ENV" >&2
    exit 1
}

[[ -d "$ASSEMBLY_DIR" ]] || {
    echo "ERROR: Assembly directory not found: $ASSEMBLY_DIR" >&2
    exit 1
}

[[ -d "$REFERENCE_DIR" ]] || {
    echo "ERROR: Reference directory not found: $REFERENCE_DIR" >&2
    exit 1
}

[[ -s "$MAPPING_FILE" ]] || {
    echo "ERROR: Mapping file missing or empty: $MAPPING_FILE" >&2
    exit 1
}

LOG="$OUTPUT_DIR/alignment.log"
exec > >(tee -a "$LOG") 2>&1

echo "=================================================================="
echo " Assembly-to-reference alignment"
echo "=================================================================="
echo "PROJECT_ROOT : $PROJECT_ROOT"
echo "Assemblies   : $ASSEMBLY_DIR"
echo "References   : $REFERENCE_DIR"
echo "Mapping file : $MAPPING_FILE"
echo "Output       : $OUTPUT_DIR"
echo "Conda env    : $CONDA_ENV"
echo "Threads      : $THREADS"
echo "Started      : $(date)"
echo "=================================================================="
echo

# =====================================================================
# Load SRR -> reference assembly mappings
# Expected columns: SRR ID, primary assembly accession, secondary accession
# =====================================================================
declare -A SRR_TO_ASSEMBLY
declare -A SRR_TO_ASSEMBLY2

while IFS=$'\t' read -r srr asm1 asm2 _rest; do
    srr="${srr//$'\r'/}"
    asm1="${asm1//$'\r'/}"
    asm2="${asm2//$'\r'/}"

    [[ -z "$srr" ]] && continue
    [[ "$srr" == "SRR ID" || "$srr" == "SRR_ID" ]] && continue

    SRR_TO_ASSEMBLY["$srr"]="$asm1"
    SRR_TO_ASSEMBLY2["$srr"]="$asm2"
done < "$MAPPING_FILE"

echo "Loaded ${#SRR_TO_ASSEMBLY[@]} SRR reference mappings"
echo

# Best_Assemblies uses flat platform_species directories:
#   ONT_ecoli
#   ONT_kpneumoniae
#   PacBio_ecoli
#   PacBio_kpneumoniae
PLATFORMS=("ONT" "PacBio")
SPECIES_LIST=("ecoli")

for PLATFORM in "${PLATFORMS[@]}"; do
    for SPECIES in "${SPECIES_LIST[@]}"; do
        ASM_BASE="$ASSEMBLY_DIR/$PLATFORM/$SPECIES"

        # ReferenceFASTA uses pneumoniae rather than kpneumoniae.
        case "$SPECIES" in
            ecoli)
                REF_SPECIES="ecoli"
                ;;
            kpneumoniae)
                REF_SPECIES="kpneumoniae"
                ;;
            *)
                echo "WARNING: Unsupported species: $SPECIES"
                continue
                ;;
        esac

        if [[ ! -d "$ASM_BASE" ]]; then
            echo "WARNING: Assembly directory not found, skipping: $ASM_BASE"
            continue
        fi

        REF_BASE="$REFERENCE_DIR/$REF_SPECIES"
        if [[ ! -d "$REF_BASE" ]]; then
            echo "WARNING: Reference species directory not found, skipping: $REF_BASE"
            continue
        fi

        # Selected assemblies are named directly as SRR/ERR/DRR accessions,
        # for example SRR12123271.fasta.
        mapfile -t ASSEMBLIES < <(
            find "$ASM_BASE" -maxdepth 1 -type f -iname '*.fasta' | sort
        )

        if [[ ${#ASSEMBLIES[@]} -eq 0 ]]; then
            echo "WARNING: No assembly FASTA files found in: $ASM_BASE"
            continue
        fi

        OUT_DIR="$OUTPUT_DIR/$PLATFORM/$SPECIES"
        mkdir -p "$OUT_DIR"

        echo "------------------------------------------------------------------"
        echo "Platform: $PLATFORM | Species: $SPECIES"
        echo "Assembly directory : $ASM_BASE"
        echo "Reference directory: $REF_BASE"
        echo "Assemblies found   : ${#ASSEMBLIES[@]}"
        echo "------------------------------------------------------------------"

        for asm_file in "${ASSEMBLIES[@]}"; do
            filename=$(basename "$asm_file")

            # Extract an SRR/ERR/DRR accession from the filename or path.
            SRR=$(printf '%s\n' "$asm_file" | grep -oE '[SED]RR[0-9]+' | head -n 1 || true)

            if [[ -z "$SRR" ]]; then
                echo "WARNING: Could not identify run accession from: $asm_file"
                continue
            fi

            ASSEMBLY="${SRR_TO_ASSEMBLY[$SRR]:-}"
            ASSEMBLY2="${SRR_TO_ASSEMBLY2[$SRR]:-}"

            if [[ -z "$ASSEMBLY" ]]; then
                echo "WARNING: No reference mapping for $SRR; skipping $filename"
                continue
            fi

            # Locate the full reference genome using either mapped accession.
            REF_GENOME=$(find "$REF_BASE" -type f \
                -path "*/$ASSEMBLY/${ASSEMBLY}*_genomic.fna" \
                -print -quit 2>/dev/null || true)

            if [[ -z "$REF_GENOME" && -n "$ASSEMBLY2" ]]; then
                REF_GENOME=$(find "$REF_BASE" -type f \
                    -path "*/$ASSEMBLY2/${ASSEMBLY2}*_genomic.fna" \
                    -print -quit 2>/dev/null || true)
            fi

            if [[ -z "$REF_GENOME" ]]; then
                echo "WARNING: Reference genome not found for $SRR ($ASSEMBLY / $ASSEMBLY2)"
                continue
            fi

            OUT_PAF="$OUT_DIR/${SRR}_assembly_vs_ref.paf"
            OUT_PAF_SR="$OUT_DIR/${SRR}_assembly_vs_ref.sr.paf"

            if [[ -s "$OUT_PAF" ]]; then
                echo "[$SRR] Standard alignment already exists; skipping"
            else
                echo "[$SRR] Standard alignment (-x asm10)"
                minimap2 \
                    -x asm10 \
                    --secondary=no \
                    -t "$THREADS" \
                    "$REF_GENOME" \
                    "$asm_file" \
                    > "$OUT_PAF"

                echo "[$SRR] Standard PAF: $(wc -l < "$OUT_PAF") alignments"
            fi

            if [[ -s "$OUT_PAF_SR" ]]; then
                echo "[$SRR] High-sensitivity alignment already exists; skipping"
            else
                echo "[$SRR] High-sensitivity alignment (-x sr)"
                minimap2 \
                    -x sr \
                    -A2 -B4 -O4,24 -E2,1 \
                    --secondary=no \
                    -t "$THREADS" \
                    "$REF_GENOME" \
                    "$asm_file" \
                    > "$OUT_PAF_SR"

                echo "[$SRR] Sensitive PAF: $(wc -l < "$OUT_PAF_SR") alignments"
            fi
        done
    done
done

echo
echo "=== Alignment completed: $(date) ==="
