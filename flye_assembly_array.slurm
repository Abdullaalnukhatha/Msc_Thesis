#!/usr/bin/env bash
#SBATCH --job-name=flye_array
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --time=48:00:00
#SBATCH --output=/rds/projects/e/elhamsak-mbru-amr/tempdata/flye_%A_%a.out
#SBATCH --error=/rds/projects/e/elhamsak-mbru-amr/tempdata/flye_%A_%a.err

set -euo pipefail

PROJ="/rds/projects/e/elhamsak-mbru-amr"
READS_DIR="$PROJ/results/ProcessedFASTQ"
OUTPUT_DIR="$PROJ/results/assemblies"
MANIFEST="$PROJ/tempdata/flye_array_manifest.tsv"
MINIFORGE="$HOME/miniforge3"
FLYE_ENV="$MINIFORGE/envs/flye_env"
THREADS="${SLURM_CPUS_PER_TASK:-12}"
ONT_Q_THRESHOLD=10

mkdir -p "$PROJ/tempdata" "$OUTPUT_DIR"

# Run this file normally to prepare the manifest and submit the array:
#   bash flye_assembly_array.slurm
# Inside an array task, SLURM_ARRAY_TASK_ID is already defined.
if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    : > "$MANIFEST"

    for platform in ONT PacBio; do
        platform_dir="$READS_DIR/$platform"
        [[ -d "$platform_dir" ]] || continue

        while IFS= read -r fastq; do
            species=$(basename "$(dirname "$fastq")")
            printf '%s\t%s\t%s\n' "$fastq" "$platform" "$species" >> "$MANIFEST"
        done < <(find "$platform_dir" -type f -name '*.fastq.gz' | sort)
    done

    total=$(wc -l < "$MANIFEST")
    if (( total == 0 )); then
        echo "ERROR: No processed FASTQ files found under: $READS_DIR" >&2
        exit 1
    fi

    echo "Created manifest with $total samples: $MANIFEST"
    echo "Submitting Flye array: 1-$total (maximum 20 concurrent tasks)"
    sbatch --array="1-${total}%20" "$0"
    exit 0
fi

module purge
module load bluebear
module load bear-apps/2024a/live
module load Miniforge3/25.3.0-3

source "$MINIFORGE/etc/profile.d/conda.sh"

if [[ -d "$FLYE_ENV" ]]; then
    conda activate "$FLYE_ENV"
elif conda env list | awk '{print $1}' | grep -qx 'flye_env'; then
    conda activate flye_env
else
    echo "ERROR: Flye environment not found at $FLYE_ENV and no named flye_env exists." >&2
    echo "Available environments:" >&2
    conda env list >&2
    exit 1
fi

for command in flye seqkit gzip awk sort uniq; do
    command -v "$command" >/dev/null 2>&1 || {
        echo "ERROR: Required command not found in the active environment: $command" >&2
        exit 1
    }
done

if [[ ! -s "$MANIFEST" ]]; then
    echo "ERROR: Manifest is missing or empty: $MANIFEST" >&2
    exit 1
fi

line=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$MANIFEST")
if [[ -z "$line" ]]; then
    echo "ERROR: No manifest entry for array task $SLURM_ARRAY_TASK_ID" >&2
    exit 1
fi

IFS=$'\t' read -r fastq platform species <<< "$line"

if [[ ! -f "$fastq" ]]; then
    echo "ERROR: FASTQ not found: $fastq" >&2
    exit 1
fi

case "$species" in
    ecoli) genome_size="5m" ;;
    kpneumoniae|pneumoniae) genome_size="5.5m" ;;
    *)
        echo "ERROR: Unknown species directory '$species' for $fastq" >&2
        exit 1
        ;;
esac

base=$(basename "$fastq" .fastq.gz)
sample=${base%%_*}
final_dir="$OUTPUT_DIR/$platform/$species"
final_fasta="$final_dir/${sample}_contigs.fasta"
sample_out="$final_dir/.flye_work_${sample}_${SLURM_ARRAY_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
dedup_fastq="$sample_out/${sample}_dedup.fastq.gz"
renamed_fastq="$sample_out/${sample}_renamed.fastq.gz"

mkdir -p "$final_dir"

if [[ -s "$final_fasta" ]]; then
    echo "Assembly already exists; skipping: $final_fasta"
    exit 0
fi

cleanup() {
    rm -rf "$sample_out"
}
trap cleanup EXIT

if ! gzip -t "$fastq" 2>/dev/null; then
    echo "ERROR: Corrupted gzip file: $fastq" >&2
    exit 1
fi

get_seqkit_value() {
    local table="$1"
    local column="$2"

    awk -F'\t' -v col="$column" '
    NR == 1 {
        for (i = 1; i <= NF; i++) {
            gsub(/^[[:space:]]+|[[:space:]]+$/, "", $i)
            if ($i == col) { idx = i; break }
        }
    }
    NR == 2 && idx {
        gsub(/,/, "", $idx)
        gsub(/^[[:space:]]+|[[:space:]]+$/, "", $idx)
        print $idx
    }' <<< "$table"
}

stats_tsv=$(seqkit stats -a -T "$fastq")
mean_len=$(get_seqkit_value "$stats_tsv" "avg_len")
mean_qual=$(get_seqkit_value "$stats_tsv" "AvgQual")
mean_len=${mean_len:-0}
mean_qual=${mean_qual:-0}

[[ "$mean_len" =~ ^[0-9]+([.][0-9]+)?$ ]] || mean_len=0
[[ "$mean_qual" =~ ^[0-9]+([.][0-9]+)?$ ]] || mean_qual=0

case "$platform" in
    ONT)
        if awk -v q="$mean_qual" -v t="$ONT_Q_THRESHOLD" 'BEGIN {exit !(q >= t)}'; then
            preset="--nano-hq"
        else
            preset="--nano-raw"
        fi
        ;;
    PacBio)
        if awk -v q="$mean_qual" 'BEGIN {exit !(q >= 20)}'; then
            preset="--pacbio-hifi"
        else
            preset="--pacbio-raw"
        fi
        ;;
    *)
        echo "ERROR: Unknown platform: $platform" >&2
        exit 1
        ;;
esac

echo "============================================================"
echo "Array job/task : ${SLURM_ARRAY_JOB_ID}/${SLURM_ARRAY_TASK_ID}"
echo "Sample         : $sample"
echo "Platform       : $platform"
echo "Species        : $species"
echo "FASTQ          : $fastq"
echo "Mean length    : $mean_len bp"
echo "Mean quality   : Q$mean_qual"
echo "Flye preset    : $preset"
echo "Genome size    : $genome_size"
echo "Threads        : $THREADS"
echo "============================================================"

mkdir -p "$sample_out"

seqkit rmdup -n "$fastq" \
    -o "$dedup_fastq" \
    --threads "$THREADS"

assembly_fastq="$dedup_fastq"
dup_check=$(gzip -cd "$dedup_fastq" | awk 'NR % 4 == 1 {print $1}' | sort | uniq -d | head -n 1 || true)

if [[ -n "$dup_check" ]]; then
    echo "Duplicate primary IDs remain; renaming headers sequentially."
    seqkit replace -p '.*' -r 'read_{nr}' "$dedup_fastq" \
        -o "$renamed_fastq" \
        --threads "$THREADS"
    assembly_fastq="$renamed_fastq"
fi

long_reads=$(seqkit seq -m 1000 -g "$assembly_fastq" 2>/dev/null \
    | seqkit stats -T \
    | awk -F'\t' 'NR==1 {for(i=1;i<=NF;i++) if($i=="num_seqs") c=i} NR==2 && c {gsub(/,/,"",$c); print $c}')
long_reads=${long_reads:-0}

if ! [[ "$long_reads" =~ ^[0-9]+$ ]] || (( long_reads == 0 )); then
    echo "ERROR: No reads of at least 1000 bp remain for $sample" >&2
    exit 1
fi

echo "Reads >=1000 bp: $long_reads"
echo "Starting Flye..."

flye "$preset" "$assembly_fastq" \
    --out-dir "$sample_out" \
    --threads "$THREADS" \
    --genome-size "$genome_size"

if [[ ! -s "$sample_out/assembly.fasta" ]]; then
    echo "ERROR: Flye finished without a non-empty assembly.fasta" >&2
    exit 1
fi

# Copy to temporary final names, then rename atomically to avoid partial outputs.
cp "$sample_out/assembly.fasta" "$final_fasta.tmp"
mv -f "$final_fasta.tmp" "$final_fasta"

if [[ -f "$sample_out/assembly_info.txt" ]]; then
    cp "$sample_out/assembly_info.txt" "$final_dir/${sample}_assembly_info.txt.tmp"
    mv -f "$final_dir/${sample}_assembly_info.txt.tmp" "$final_dir/${sample}_assembly_info.txt"
fi

if [[ -f "$sample_out/assembly_graph.gfa" ]]; then
    cp "$sample_out/assembly_graph.gfa" "$final_dir/${sample}_assembly_graph.gfa.tmp"
    mv -f "$final_dir/${sample}_assembly_graph.gfa.tmp" "$final_dir/${sample}_assembly_graph.gfa"
fi

echo "Completed: $final_fasta"
