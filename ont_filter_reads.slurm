#!/bin/bash
#SBATCH --job-name=ONT_filter
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --output=/rds/projects/e/elhamsak-mbru-amr/tempdata/ONT_filter_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-mbru-amr/tempdata/ONT_filter_%j.err

set -euo pipefail

module purge
module load bluebear
module load bear-apps/2024a/live
module load Miniforge3/25.3.0-3

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate /rds/homes/a/axa2273/miniforge3/envs/porechop_env

PROJ="/rds/projects/e/elhamsak-mbru-amr"
INPUT_DIR="${1:-$PROJ/data/FASTQ}"
OUTPUT_DIR="${2:-$PROJ/results/ProcessedFASTQ/ONT}"
MIN_LENGTH="${3:-500}"
MIN_QUALITY="${4:-10}"
THREADS="${SLURM_CPUS_PER_TASK:-12}"

mkdir -p "$OUTPUT_DIR" "$PROJ/tempdata"

for cmd in porechop chopper seqkit gzip bc; do
    command -v "$cmd" >/dev/null 2>&1 || {
        echo "ERROR: Required command not found: $cmd" >&2
        exit 1
    }
done

echo "=================================================================="
echo " ONT read filtering"
echo "=================================================================="
echo "Job ID       : ${SLURM_JOB_ID:-N/A}"
echo "Node         : $(hostname)"
echo "Input        : $INPUT_DIR"
echo "Output       : $OUTPUT_DIR"
echo "Min length   : $MIN_LENGTH bp"
echo "Min quality  : Q$MIN_QUALITY"
echo "Threads      : $THREADS"
echo "Started      : $(date)"
echo "=================================================================="

SPECIES_LIST=("ecoli" "kpneumoniae")

for SPECIES in "${SPECIES_LIST[@]}"; do
    SPECIES_INPUT="$INPUT_DIR/$SPECIES"
    SPECIES_OUTPUT="$OUTPUT_DIR/$SPECIES"

    echo "=== Processing species: $SPECIES ==="

    if [[ ! -d "$SPECIES_INPUT" ]]; then
        echo "WARNING: Directory not found, skipping: $SPECIES_INPUT"
        continue
    fi

    mkdir -p "$SPECIES_OUTPUT"

    # Restrict the search to ONT directories so PacBio files are not processed.
    mapfile -t FASTQ_FILES < <(
        find "$SPECIES_INPUT" -type f \
            \( -path "*/ONT/*.fastq.gz" -o -path "*/ONT/*.fq.gz" \) | sort
    )

    TOTAL_FILES=${#FASTQ_FILES[@]}
    if (( TOTAL_FILES == 0 )); then
        echo "WARNING: No ONT FASTQ files found under $SPECIES_INPUT, skipping"
        continue
    fi

    echo "Found $TOTAL_FILES ONT FASTQ files"

    STATS_FILE="$SPECIES_OUTPUT/${SPECIES}_processing_stats.tsv"
    echo -e "Sample\tInput_Reads\tInput_Bases\tFiltered_Reads\tFiltered_Bases\tRetention_Rate" > "$STATS_FILE"

    CURRENT=0
    for fastq in "${FASTQ_FILES[@]}"; do
        ((CURRENT+=1))
        SAMPLE=$(basename "$fastq" | sed -E 's/\.(fastq|fq)\.gz$//')
        OUT_FASTQ="$SPECIES_OUTPUT/${SAMPLE}_processed.fastq.gz"
        TMP_FASTQ="${OUT_FASTQ}.tmp"

        echo "[$CURRENT/$TOTAL_FILES] Processing: $SAMPLE"

        if ! gzip -t "$fastq" 2>/dev/null; then
            echo "WARNING: Corrupted file, skipping: $fastq"
            continue
        fi

        INPUT_STATS=$(seqkit stats "$fastq" -T | tail -n1)
        INPUT_READS=$(echo "$INPUT_STATS" | cut -f4 | tr -d ',')
        INPUT_BASES=$(echo "$INPUT_STATS" | cut -f5 | tr -d ',')

        rm -f "$TMP_FASTQ"
        porechop -i "$fastq" --threads "$THREADS" --discard_middle --verbosity 0 \
            | chopper --quality "$MIN_QUALITY" --minlength "$MIN_LENGTH" --threads "$THREADS" \
            | gzip -c > "$TMP_FASTQ"

        gzip -t "$TMP_FASTQ"
        mv -f "$TMP_FASTQ" "$OUT_FASTQ"

        OUTPUT_STATS=$(seqkit stats "$OUT_FASTQ" -T | tail -n1)
        OUTPUT_READS=$(echo "$OUTPUT_STATS" | cut -f4 | tr -d ',')
        OUTPUT_BASES=$(echo "$OUTPUT_STATS" | cut -f5 | tr -d ',')

        if (( INPUT_READS > 0 )); then
            RETENTION=$(echo "scale=2; $OUTPUT_READS * 100 / $INPUT_READS" | bc)
        else
            RETENTION="0.00"
        fi

        echo -e "$SAMPLE\t$INPUT_READS\t$INPUT_BASES\t$OUTPUT_READS\t$OUTPUT_BASES\t${RETENTION}%" >> "$STATS_FILE"
        echo "Retained: $OUTPUT_READS / $INPUT_READS reads (${RETENTION}%)"
    done

done

echo "=================================================================="
echo "ONT filtering completed: $(date)"
echo "=================================================================="
