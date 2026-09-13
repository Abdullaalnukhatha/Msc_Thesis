#!/bin/bash
#SBATCH --job-name=ONT_QC
#SBATCH --time=24:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --output=/rds/projects/e/elhamsak-mbru-amr/tempdata/ont_qc_%j.out

set -euo pipefail
export LC_ALL=C

module purge
module load bluebear
module load bear-apps/2024a/live
module load Miniforge3/25.3.0-3

source /rds/homes/a/axa2273/miniforge3/etc/profile.d/conda.sh

INPUT_DIR="/rds/projects/e/elhamsak-mbru-amr/data/FASTQ"
OUTPUT_DIR="/rds/projects/e/elhamsak-mbru-amr/results/qc/ONT"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

mkdir -p "$OUTPUT_DIR"

echo "=== ONT Quality Control Pipeline ==="
echo "Input directory: $INPUT_DIR"
echo "Output directory: $OUTPUT_DIR"
echo "Threads: $THREADS"
echo ""

SPECIES_LIST=("ecoli" "kpneumoniae")

for SPECIES in "${SPECIES_LIST[@]}"; do
    SPECIES_INPUT="$INPUT_DIR/$SPECIES"
    SPECIES_OUTPUT="$OUTPUT_DIR/$SPECIES"

    echo "=== Processing species: $SPECIES ==="

    if [ ! -d "$SPECIES_INPUT" ]; then
        echo "WARNING: Directory not found, skipping: $SPECIES_INPUT"
        continue
    fi

    mkdir -p "$SPECIES_OUTPUT"

    mapfile -t FASTQ_FILES < <(find "$SPECIES_INPUT" -type f -path "*/ONT/*.fastq.gz" | sort)
    TOTAL_FILES="${#FASTQ_FILES[@]}"

    if [ "$TOTAL_FILES" -eq 0 ]; then
        echo "WARNING: No FASTQ files found under $SPECIES_INPUT/*/ONT/*.fastq.gz, skipping"
        continue
    fi

    echo "Found $TOTAL_FILES FASTQ files"
    echo ""

    conda activate nanoplot_env

    CURRENT=0
    for fastq in "${FASTQ_FILES[@]}"; do
        CURRENT=$((CURRENT + 1))
        SAMPLE=$(basename "$fastq")
        SAMPLE="${SAMPLE%.fastq.gz}"
        SAMPLE="${SAMPLE%.fq.gz}"

        SAMPLE_OUT="$SPECIES_OUTPUT/$SAMPLE"
        REPORT_HTML="$SAMPLE_OUT/${SAMPLE}_NanoPlot-report.html"

        if [ -f "$REPORT_HTML" ]; then
            echo "[$CURRENT/$TOTAL_FILES] [SKIP] NanoPlot already done for $SAMPLE"
            continue
        fi

        echo "[$CURRENT/$TOTAL_FILES] Processing: $SAMPLE"

        if ! gzip -t "$fastq" 2>/dev/null; then
            echo "WARNING: Corrupted file, skipping: $fastq"
            continue
        fi

        mkdir -p "$SAMPLE_OUT"

        NanoPlot \
            --fastq "$fastq" \
            --outdir "$SAMPLE_OUT" \
            --prefix "${SAMPLE}_" \
            --threads "$THREADS" \
            --N50 \
            --only-report

        echo ""
    done

    conda activate multiqc_env

    MULTIQC_REPORT="$SPECIES_OUTPUT/${SPECIES}_ont_multiqc_report.html"
    if [ -f "$MULTIQC_REPORT" ]; then
        echo "[SKIP] MultiQC already exists for $SPECIES"
    else
        echo "Generating MultiQC summary..."
        multiqc "$SPECIES_OUTPUT" -o "$SPECIES_OUTPUT" -n "${SPECIES}_ont_multiqc_report" --force --quiet
    fi

    conda activate nanoplot_env

    SUMMARY_FILE="$SPECIES_OUTPUT/${SPECIES}_qc_summary.tsv"
    echo "Creating summary table..."
    echo -e "Sample\tTotal_Reads\tTotal_Bases\tMean_Length\tMedian_Length\tN50\tMean_Quality\tMedian_Quality\tCoverage_5.5Mb" > "$SUMMARY_FILE"

    find "$SPECIES_OUTPUT" -type f \( -name "*_NanoStats.txt" -o -name "*_NanoPlot_stats.txt" \) | sort | while read -r stats_file; do
        SAMPLE=$(basename "$(dirname "$stats_file")")

        TOTAL_READS=$(grep -m1 "Number of reads:" "$stats_file" | awk '{print $NF}' | tr -d ',' || true)
        TOTAL_BASES=$(grep -m1 "Total bases:" "$stats_file" | awk '{print $NF}' | tr -d ',' || true)
        MEAN_LENGTH=$(grep -m1 "Mean read length:" "$stats_file" | awk '{print $NF}' | tr -d ',' || true)
        MEDIAN_LENGTH=$(grep -m1 "Median read length:" "$stats_file" | awk '{print $NF}' | tr -d ',' || true)
        N50=$(grep -m1 "Read length N50:" "$stats_file" | awk '{print $NF}' | tr -d ',' || true)
        MEAN_QUAL=$(grep -m1 "Mean read quality:" "$stats_file" | awk '{print $NF}' || true)
        MEDIAN_QUAL=$(grep -m1 "Median read quality:" "$stats_file" | awk '{print $NF}' || true)

        TOTAL_BASES="${TOTAL_BASES:-0}"
        COVERAGE=$(echo "scale=1; $TOTAL_BASES / 5500000" | bc)

        echo -e "$SAMPLE\t${TOTAL_READS:-NA}\t${TOTAL_BASES:-NA}\t${MEAN_LENGTH:-NA}\t${MEDIAN_LENGTH:-NA}\t${N50:-NA}\t${MEAN_QUAL:-NA}\t${MEDIAN_QUAL:-NA}\t${COVERAGE}x"
    done >> "$SUMMARY_FILE"

    echo ""
    echo "=== Quality Control Complete: $SPECIES ==="
    echo "MultiQC report: $SPECIES_OUTPUT/${SPECIES}_ont_multiqc_report.html"
    echo "Summary table:  $SUMMARY_FILE"
    echo ""
done

echo "ONT QC finished at $(date)"
