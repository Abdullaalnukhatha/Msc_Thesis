#!/bin/bash
#SBATCH --job-name=PacBio_QC
#SBATCH --time=12:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=4
#SBATCH --output=/rds/projects/e/elhamsak-mbru-amr/tempdata/pacbio_qc_%j.out

set -euo pipefail
export LC_ALL=C

module purge
module load bluebear
module load bear-apps/2024a/live
module load Miniforge3/25.3.0-3

source /rds/homes/a/axa2273/miniforge3/etc/profile.d/conda.sh
conda activate porechop_env

INPUT_DIR="/rds/projects/e/elhamsak-mbru-amr/data/FASTQ"
OUTPUT_DIR="/rds/projects/e/elhamsak-mbru-amr/results/qc/PacBio"
THREADS="${SLURM_CPUS_PER_TASK:-4}"

mkdir -p "$OUTPUT_DIR"

echo "=== PacBio Quality Control Pipeline ==="
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

    mapfile -t FASTQ_FILES < <(find "$SPECIES_INPUT" -type f -path "*/PacBio/*.fastq.gz" | sort)
    TOTAL_FILES="${#FASTQ_FILES[@]}"

    if [ "$TOTAL_FILES" -eq 0 ]; then
        echo "WARNING: No FASTQ files found under $SPECIES_INPUT/*/PacBio/*.fastq.gz, skipping"
        continue
    fi

    echo "Found $TOTAL_FILES FASTQ files"
    echo ""

    CURRENT=0
    for fastq in "${FASTQ_FILES[@]}"; do
        CURRENT=$((CURRENT + 1))
        SAMPLE=$(basename "$fastq")
        SAMPLE="${SAMPLE%.fastq.gz}"
        SAMPLE="${SAMPLE%.fq.gz}"

        STATS_OUT="$SPECIES_OUTPUT/${SAMPLE}_seqkit_stats.txt"

        if [ -f "$STATS_OUT" ]; then
            echo "[$CURRENT/$TOTAL_FILES] [SKIP] PacBio QC already done for $SAMPLE"
            continue
        fi

        echo "[$CURRENT/$TOTAL_FILES] Processing: $SAMPLE"

        if ! gzip -t "$fastq" 2>/dev/null; then
            echo "WARNING: Corrupted file, skipping: $fastq"
            continue
        fi

        seqkit stats "$fastq" -T -j "$THREADS" > "$STATS_OUT"
        echo "  Completed"
    done

    echo "Creating summary table..."

    SUMMARY_TSV="$SPECIES_OUTPUT/${SPECIES}_qc_summary.tsv"

    if [ "$SPECIES" = "ecoli" ]; then
        GENOME_SIZE=5000000
    else
        GENOME_SIZE=5500000
    fi

    echo -e "sample\tfile\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len\tcoverage" > "$SUMMARY_TSV"

    mapfile -t STATS_FILES < <(find "$SPECIES_OUTPUT" -maxdepth 1 -name "*_seqkit_stats.txt" | sort)

    if [ "${#STATS_FILES[@]}" -eq 0 ]; then
        echo "WARNING: No seqkit stats files found for $SPECIES, skipping summary table population"
    else
        for stats_file in "${STATS_FILES[@]}"; do
            SAMPLE=$(basename "$stats_file" _seqkit_stats.txt)

            row=$(tail -n +2 "$stats_file" | head -n 1 || true)
            [ -n "$row" ] || continue

            file=$(echo "$row" | cut -f1)
            format=$(echo "$row" | cut -f2)
            type=$(echo "$row" | cut -f3)
            num_seqs=$(echo "$row" | cut -f4 | tr -d ',')
            sum_len=$(echo "$row" | cut -f5 | tr -d ',')
            min_len=$(echo "$row" | cut -f6 | tr -d ',')
            avg_len=$(echo "$row" | cut -f7 | tr -d ',')
            max_len=$(echo "$row" | cut -f8 | tr -d ',')

            coverage=$(echo "scale=2; ${sum_len:-0} / $GENOME_SIZE" | bc)

            echo -e "${SAMPLE}\t${file}\t${format}\t${type}\t${num_seqs}\t${sum_len}\t${min_len}\t${avg_len}\t${max_len}\t${coverage}"
        done >> "$SUMMARY_TSV"
    fi

    echo ""
    echo "=== Quality Control Complete: $SPECIES ==="
    echo "Summary table: $SUMMARY_TSV"
    echo "Genome size used for coverage: $GENOME_SIZE"
    echo ""
done

echo "PacBio QC finished at $(date)"
