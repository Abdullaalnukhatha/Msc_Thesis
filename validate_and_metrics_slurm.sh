#!/bin/bash
# validate_and_metrics.sh

# --- SLURM Batch Directives (Ignored when run locally) -----------------
#SBATCH --job-name=validate_metrics
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=24:00:00
#SBATCH --output=validate_metrics_%j.out
#SBATCH --error=validate_metrics_%j.err

set -Eeuo pipefail

#module purge
#module load bluebear
#module load bear-apps/2024a/live

trap 'echo "ERROR at line $LINENO: $BASH_COMMAND" >&2' ERR

# --- Define PROJECT_ROOT based on environment -------------------------
if [ -n "${SLURM_JOB_ID:-}" ]; then
    # Cluster/SLURM execution environment
    PROJECT_ROOT="${PROJECT_ROOT:-/rds/projects/e/elhamsak-mbru-amr}"
    CONDA_ENV="${CONDA_ENV:-$PROJECT_ROOT/envs/ref_valid}"

    # --- Arguments --------------------------------------------------------
    PREDICTIONS_DIR="${1:-$PROJECT_ROOT/results/plasmid_benchmark}"
    REFERENCE_DIR="${2:-$PROJECT_ROOT/data/ReferenceFASTA}"
    ASSEMBLY_DIR="${3:-$PROJECT_ROOT/results/Best_Assemblies}"
    ALIGNMENT_DIR="${4:-$PROJECT_ROOT/results/alignments}"
    MAPPING_FILE="${5:-$PROJECT_ROOT/data/ReferenceFASTA/reference_mapping.tsv}"
    OUTPUT_DIR="${6:-$PROJECT_ROOT/results/plasmid_benchmark/validation_metrics}"
    SCRIPTS_DIR="${7:-$PROJECT_ROOT/Accepted_scripts}"
    THREADS="${8:-${SLURM_CPUS_PER_TASK:-8}}"
    MIN_IDENTITY="${9:-95}"

    SPECIES_LIST=("ecoli")
        
    # Activate Conda environment if path exists
    if command -v conda >/dev/null 2>&1 && [ -d "$CONDA_ENV" ]; then
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate "$CONDA_ENV"
    fi
else
    # Local execution environment (Current Working Directory)
    PROJECT_ROOT="."
    PREDICTIONS_DIR="${1:-$PROJECT_ROOT/results/plasmid_benchmark}"
    REFERENCE_DIR="${2:-$PROJECT_ROOT/data/FASTA}"
    ASSEMBLY_DIR="${3:-$PROJECT_ROOT/results/assemblies}"
    ALIGNMENT_DIR="${4:-$PROJECT_ROOT/results/alignments}"
    MAPPING_FILE="${5:-$PROJECT_ROOT/data/reference_mapping.tsv}"
    OUTPUT_DIR="${6:-$PROJECT_ROOT/results/validation}"
    SCRIPTS_DIR="${7:-$PROJECT_ROOT/scripts}"
    THREADS="${8:-${SLURM_CPUS_PER_TASK:-8}}"
    MIN_IDENTITY="${9:-95}"

    SPECIES_LIST=("ecoli")
fi

mkdir -p "$OUTPUT_DIR"
LOG="$OUTPUT_DIR/validation.log"
exec > >(tee "$LOG") 2>&1

echo "=================================================================="
echo " Plasmid Validation & Benchmarking Metrics"
echo "=================================================================="
echo "PROJECT_ROOT : $PROJECT_ROOT"
echo "Predictions  : $PREDICTIONS_DIR"
echo "References   : $REFERENCE_DIR"
echo "Alignments   : $ALIGNMENT_DIR"
echo "Mapping file : $MAPPING_FILE"
echo "Output       : $OUTPUT_DIR"
echo "Scripts      : $SCRIPTS_DIR"
echo "Threads      : $THREADS"
echo "Min Identity : ${MIN_IDENTITY}%"
echo "Started      : $(date)"
echo "=================================================================="
echo ""

# --- Validate required files ---------------------------------------------------
for f in \
    "$MAPPING_FILE" \
    "$SCRIPTS_DIR/parse_tool_report.py" \
    "$SCRIPTS_DIR/calc_metrics.py" \
    "$SCRIPTS_DIR/get_ref_size_category.py"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required file not found: $f"
        exit 1
    fi
done

# --- Load SRR → Assembly mapping ------------------------------------------------
declare -A SRR_TO_ASSEMBLY
declare -A SRR_TO_ASSEMBLY2
while IFS=$'\t' read -r srr asm1 asm2; do
    [[ "$srr" == "SRR ID" || -z "$srr" ]] && continue
    SRR_TO_ASSEMBLY["$srr"]="$asm1"
    SRR_TO_ASSEMBLY2["$srr"]="$asm2"
done < "$MAPPING_FILE"
echo "Loaded ${#SRR_TO_ASSEMBLY[@]} SRR → Assembly mappings"
echo ""

# --- Helper: check empty FASTA ---------------------------------------------------
is_empty_fasta() {
    local f="$1"
    [ ! -f "$f" ] && return 0
    [ ! -s "$f" ] && return 0
    grep -q '^>' "$f" && return 1 || return 0
}

# --- Tool configuration ------------------------------------------------------------
# Patterns use literal {SRR} as placeholder, resolved at runtime with bash
# parameter substitution: ${pattern//\{SRR\}/$SRR}

declare -A TOOL_REPORT_PATTERN
declare -A TOOL_FASTA_PATTERN

TOOL_REPORT_PATTERN["MOB-suite"]="contig_report.txt"
TOOL_FASTA_PATTERN["MOB-suite"]="plasmid_*.fasta"

TOOL_REPORT_PATTERN["Plasmer"]="results/{SRR}.plasmer.predClass.tsv"
TOOL_FASTA_PATTERN["Plasmer"]="results/{SRR}.plasmer.predPlasmids.fa"

TOOL_REPORT_PATTERN["Plassembler"]="run/plassembler_summary.tsv"
TOOL_FASTA_PATTERN["Plassembler"]="run/plassembler_plasmids.fasta"

TOOL_REPORT_PATTERN["PlasmidEC"]="run/ensemble_output.csv"
TOOL_FASTA_PATTERN["PlasmidEC"]="run/plasmid_contigs.fasta"

TOOL_REPORT_PATTERN["Platon"]="run/{SRR}.tsv"
TOOL_FASTA_PATTERN["Platon"]="run/{SRR}.plasmid.fasta"

TOOL_REPORT_PATTERN["RFPlasmid"]="run/prediction_full.csv"
TOOL_FASTA_PATTERN["RFPlasmid"]="{SRR}_plasmids.fasta"

TOOLS=("MOB-suite" "Plasmer" "PlasmidEC" "Plassembler" "Platon" "RFPlasmid")
PLATFORMS=("ONT" "PacBio")

# --- Metrics output file ---------------------------------------------------------
METRICS_FILE="$OUTPUT_DIR/all_metrics.tsv"
echo -e "Tool\tPlatform\tSpecies\tSample\t"\
"Ref_Plasmid_Count\tPred_Plasmid_Count\t"\
"TP\tFP\tFN\t"\
"Accuracy\tCircular_TP\tCircularization_Pct\t"\
"Sensitivity\tContig_Precision\tBP_Precision\tF1_Score\t"\
"Mean_Identity\tMean_Completeness\tContamination_Pct\tSize_Category" \
    > "$METRICS_FILE"

# --- Main loop ---------------------------------------------------------------
for TOOL in "${TOOLS[@]}"; do
    for PLATFORM in "${PLATFORMS[@]}"; do
        for SPECIES in "${SPECIES_LIST[@]}"; do

            PRED_DIR="$PREDICTIONS_DIR/$TOOL/$PLATFORM/$SPECIES"
            if [ ! -d "$PRED_DIR" ]; then
                echo "WARNING: Not found, skipping: $PRED_DIR"
                continue
            fi

            mapfile -t SAMPLE_DIRS < <(find "$PRED_DIR" \
                -maxdepth 2 -type d | \
                grep -E '/[SED]RR[0-9]+$' | sort)

            if [ "${#SAMPLE_DIRS[@]}" -eq 0 ]; then
                echo "WARNING: No sample directories in $PRED_DIR"
                continue
            fi

            echo "----------------------------------------------------------------"
            echo "Tool: $TOOL | Platform: $PLATFORM | Species: $SPECIES"
            echo "Samples: ${#SAMPLE_DIRS[@]}"
            echo ""

            CURRENT=0
            TOTAL="${#SAMPLE_DIRS[@]}"

            for sample_dir in "${SAMPLE_DIRS[@]}"; do
                CURRENT=$((CURRENT + 1))
                SRR=$(basename "$sample_dir")
                ASSEMBLY="${SRR_TO_ASSEMBLY[$SRR]:-}"
                ASSEMBLY2="${SRR_TO_ASSEMBLY2[$SRR]:-}"

                if [ -z "$ASSEMBLY" ]; then
                    echo "[$CURRENT/$TOTAL] WARNING: No mapping for $SRR, skipping"
                    continue
                fi

                echo "[$CURRENT/$TOTAL] $SRR → $ASSEMBLY"

                # --- Locate reference files ------------------------------------------------
                SPECIES_DIRS=("$SPECIES")

                if [ "$SPECIES" = "kpneumoniae" ]; then
                    SPECIES_DIRS=("kpneumoniae" "pneumoniae")
                fi

                REF_GENOME=""
                REF_PLASMID=""

                for SPECIES_DIR in "${SPECIES_DIRS[@]}"; do
                    REF_GENOME=$(find "$REFERENCE_DIR/$SPECIES_DIR" \
                        -path "*/$ASSEMBLY/${ASSEMBLY}*_genomic.fna" \
                        2>/dev/null | head -1 || true)

                    REF_PLASMID=$(find "$REFERENCE_DIR/$SPECIES_DIR" \
                        -path "*/$ASSEMBLY/${ASSEMBLY}*_plasmid.fna" \
                        2>/dev/null | head -1 || true)

                    [ -n "$REF_GENOME" ] && break
                done

                # Fallback to secondary assembly accession
                if [ -z "$REF_GENOME" ] && [ -n "$ASSEMBLY2" ]; then
                    for SPECIES_DIR in "${SPECIES_DIRS[@]}"; do
                        REF_GENOME=$(find "$REFERENCE_DIR/$SPECIES_DIR" \
                            -path "*/$ASSEMBLY2/${ASSEMBLY2}*_genomic.fna" \
                            2>/dev/null | head -1 || true)

                        REF_PLASMID=$(find "$REFERENCE_DIR/$SPECIES_DIR" \
                            -path "*/$ASSEMBLY2/${ASSEMBLY2}*_plasmid.fna" \
                            2>/dev/null | head -1 || true)

                        [ -n "$REF_GENOME" ] && break
                    done
                fi

                if [ -z "$REF_GENOME" ]; then
                    echo "  WARNING: Reference genome not found for $ASSEMBLY, skipping"
                    continue
                fi

                # --- Check reference plasmid emptiness -------------------------------------
                REF_EMPTY=false
                { [ -z "$REF_PLASMID" ] || is_empty_fasta "$REF_PLASMID"; } \
                    && REF_EMPTY=true

                REF_COUNT=0
                [ "$REF_EMPTY" = false ] && \
                    REF_COUNT=$(grep -c '^>' "$REF_PLASMID" 2>/dev/null || echo 0)

                # --- Get size category from reference plasmid ------------------------------
                if [ "$REF_EMPTY" = false ]; then
                    REF_SIZE_CAT=$(python3 "$SCRIPTS_DIR/get_ref_size_category.py" \
                        "$REF_PLASMID")
                else
                    REF_SIZE_CAT="None"
                fi

                SAMPLE_OUT="$OUTPUT_DIR/$TOOL/$PLATFORM/$SPECIES"
                mkdir -p "$SAMPLE_OUT"

                # --- Locate pre-computed PAF ------------------------------------------------

                if [ "$TOOL" != "Plassembler" ]; then
                    PAF_FILE="$ALIGNMENT_DIR/$PLATFORM/$SPECIES/${SRR}_assembly_vs_ref.paf"
                else
                    PAF_FILE="$ALIGNMENT_DIR/$PLATFORM/$SPECIES/${SRR}_assembly_vs_ref.sr.paf"
                fi

                if [ ! -f "$PAF_FILE" ]; then
                    echo "  WARNING: PAF not found: $PAF_FILE, skipping"
                    continue
                fi

                # --- Resolve tool report and plasmid FASTA paths ----------------------------
                REPORT_PATTERN="${TOOL_REPORT_PATTERN[$TOOL]}"
                FASTA_PATTERN="${TOOL_FASTA_PATTERN[$TOOL]}"

                # Substitute {SRR} placeholder at runtime
                REPORT_FILE="${sample_dir}/${REPORT_PATTERN//\{SRR\}/$SRR}"
                PLASMID_FASTA="${sample_dir}/${FASTA_PATTERN//\{SRR\}/$SRR}"

                # Absolute path normalization guard
                if [ -n "$PLASMID_FASTA" ]; then
                    PLASMID_FASTA=$(readlink -f "$PLASMID_FASTA" 2>/dev/null || echo "$PLASMID_FASTA")
                fi

                echo " Plasmid FASTA: $PLASMID_FASTA"
                echo " Report file  : $REPORT_FILE"

                # MOB-suite: concatenate multiple plasmid_*.fasta files
                if [ "$TOOL" = "MOB-suite" ]; then
                    CONCAT_FASTA="$SAMPLE_OUT/${SRR}_concat_plasmids.fasta"
                    if ls "$sample_dir"/plasmid_*.fasta 2>/dev/null | grep -q .; then
                        cat "$sample_dir"/plasmid_*.fasta > "$CONCAT_FASTA"
                        PLASMID_FASTA="$CONCAT_FASTA"
                    else
                        PLASMID_FASTA=""
                    fi
                fi

                # Plasmer: fallback filename search
                if [ "$TOOL" = "Plasmer" ] && [ ! -f "$PLASMID_FASTA" ]; then
                    PLASMER_FALLBACK=$(find "$sample_dir/results" \
                        -name "*.predPlasmids.fa" 2>/dev/null | head -1 || true)
                    [ -n "$PLASMER_FALLBACK" ] && PLASMID_FASTA="$PLASMER_FALLBACK"
                fi

                # --- Missing report file → log as FN --------------------------------------
                if [ ! -f "$REPORT_FILE" ]; then
                    echo "  WARNING: Report file not found: $REPORT_FILE"
                    echo -e "$TOOL\t$PLATFORM\t$SPECIES\t$SRR\t"\
"$REF_COUNT\t0\t0\t0\t$REF_COUNT\t"\
"NA\t0\t0\t"\
"0\tNA\tNA\t0\t"\
"NA\t0\t0\t$REF_SIZE_CAT" \
                        >> "$METRICS_FILE"
                    continue
                fi

                # --- Parse tool report → predicted plasmid contigs --------------------------
                PRED_TSV="$SAMPLE_OUT/${SRR}_predicted_plasmids.tsv"
                python3 "$SCRIPTS_DIR/parse_tool_report.py" \
                    "$TOOL" "$REPORT_FILE" "$PLASMID_FASTA" \
                    > "$PRED_TSV"

                PRED_COUNT=$(wc -l < "$PRED_TSV" || echo 0)

                # -------------------------------- Handle edge cases --------------------------
                # Case 1: No ref plasmids, no predictions
                # -----------------------------------------------------------------------------
                if [ "$REF_EMPTY" = true ] && [ "$PRED_COUNT" -eq 0 ]; then
                    echo "  Plasmid-free isolate, no predictions"
                    echo -e "$TOOL\t$PLATFORM\t$SPECIES\t$SRR\t"\
"0\t0\t"\
"NA\tNA\tNA\t"\
"NA\tNA\tNA\t"\
"NA\tNA\tNA\tNA\t"\
"NA\tNA\t0\t$REF_SIZE_CAT" \
                        >> "$METRICS_FILE"
                    rm -f "$PRED_TSV"
                    continue
                fi

                # ----------------------------------------------------------------------------
                # Case 2: No ref plasmids, has predictions → all FP
                # ----------------------------------------------------------------------------
                if [ "$REF_EMPTY" = true ] && [ "$PRED_COUNT" -gt 0 ]; then
                    FP_SIZE=$(awk '{sum+=$2} END{print sum+0}' "$PRED_TSV")
                    echo "  FP only: $PRED_COUNT prediction(s), no reference plasmids"
                    echo -e "$TOOL\t$PLATFORM\t$SPECIES\t$SRR\t"\
"0\t$PRED_COUNT\t"\
"0\t$FP_SIZE\t0\t"\
"NA\t0\t0\t"\
"0\t0\t0\t0\t"\
"NA\tNA\t100\t$REF_SIZE_CAT" \
                        >> "$METRICS_FILE"
                    rm -f "$PRED_TSV"
                    continue
                fi

                # -----------------------------------------------------------------------------
                # Case 3: Has ref plasmids, no predictions → all FN
                # -----------------------------------------------------------------------------
                if [ "$REF_EMPTY" = false ] && [ "$PRED_COUNT" -eq 0 ]; then
                    echo "  FN only: $REF_COUNT reference plasmid(s) missed"
                    echo -e "$TOOL\t$PLATFORM\t$SPECIES\t$SRR\t"\
"$REF_COUNT\t0\t"\
"0\t0\t$REF_COUNT\t"\
"0\t0\t0\t"\
"0\t0\t0\t0\t"\
"NA\t0\t0\t$REF_SIZE_CAT" \
                        >> "$METRICS_FILE"
                    rm -f "$PRED_TSV"
                    continue
                fi

                # ------------------------------------------------------------------------------
                # Case 4: Both present → full alignment-block metric calculation
                # ------------------------------------------------------------------------------
                echo "  Metrics: $PRED_COUNT predicted vs $REF_COUNT reference"
                CIRC_PAF="$SAMPLE_OUT/${SRR}_circ.paf"

                if [ -n "$PLASMID_FASTA" ] && [ -f "$PLASMID_FASTA" ]; then
                    minimap2 \
                        -x asm5 \
                        --secondary=no \
                        -t "$THREADS" \
                        "$PLASMID_FASTA" \
                        "$PLASMID_FASTA" \
                        > "$CIRC_PAF" 2>/dev/null

                    CIRC_LINES=$(wc -l < "$CIRC_PAF" || echo 0)
                    echo "  Circ PAF alignments: $CIRC_LINES"
                else
                    touch "$CIRC_PAF"
                    echo "  Circ PAF: no plasmid FASTA available"
                fi

                # --- Calculate all metrics via Python ------------------------------------------
                METRICS=$(python3 "$SCRIPTS_DIR/calc_metrics.py" \
                    "$PAF_FILE" \
                    "$REF_GENOME" \
                    "$PRED_TSV" \
                    "$MIN_IDENTITY" \
                    "$CIRC_PAF" \
                    "$REF_SIZE_CAT")

                echo "  RAW METRICS: $METRICS"

                # Extract fields (14 total from calc_metrics.py)
                TP=$(echo "$METRICS"              | cut -f1)
                FP=$(echo "$METRICS"              | cut -f2)
                FN=$(echo "$METRICS"              | cut -f3)
                ACCURACY=$(echo "$METRICS"        | cut -f4)
                CIRC_TP=$(echo "$METRICS"         | cut -f5)
                CIRC_PCT=$(echo "$METRICS"        | cut -f6)
                SENS=$(echo "$METRICS"            | cut -f7)
                PREC_C=$(echo "$METRICS"          | cut -f8)
                PREC_BP=$(echo "$METRICS"         | cut -f9)
                F1=$(echo "$METRICS"              | cut -f10)
                IDENTITY=$(echo "$METRICS"        | cut -f11)
                COMPLETENESS=$(echo "$METRICS"    | cut -f12)
                CONTAMINATION=$(echo "$METRICS"   | cut -f13)
                SIZE_CAT=$(echo "$METRICS"        | cut -f14)

                # Guard: if SIZE_CAT still empty, use REF_SIZE_CAT
                [ -z "$SIZE_CAT" ] && SIZE_CAT="$REF_SIZE_CAT"

                # Write metrics row
                echo -e "$TOOL\t$PLATFORM\t$SPECIES\t$SRR\t"\
"$REF_COUNT\t$PRED_COUNT\t"\
"$TP\t$FP\t$FN\t"\
"$ACCURACY\t$CIRC_TP\t$CIRC_PCT\t"\
"$SENS\t$PREC_C\t$PREC_BP\t$F1\t"\
"$IDENTITY\t$COMPLETENESS\t$CONTAMINATION\t$SIZE_CAT" \
                    >> "$METRICS_FILE"

                echo "  TP=$TP | FP=$FP | FN=$FN"
                echo "  F1=$F1 | Sensitivity=$SENS | BP_Precision=$PREC_BP"
                echo "  Completeness=${COMPLETENESS}% | Contamination=${CONTAMINATION}%"
                echo "  Circularization=${CIRC_PCT}% | Size=${SIZE_CAT}"

                # --- Cleanup intermediates ----------------------------------------------------
                rm -f "$CIRC_PAF" "$PRED_TSV"
                echo ""

            done
        done
    done
done

# --- Aggregate results -------------------------------------------------
echo "=================================================================="
echo " Aggregating Results"
echo "=================================================================="
python3 "$SCRIPTS_DIR/aggregate_metrics.py" "$METRICS_FILE"

echo ""
echo "=================================================================="
echo " Validation Complete"
echo "=================================================================="
echo "Metrics : $METRICS_FILE"
echo "Log     : $LOG"
echo "Finished: $(date)"
