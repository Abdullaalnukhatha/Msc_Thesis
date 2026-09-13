#!/bin/bash

set -euo pipefail

PROJECT="/rds/projects/e/elhamsak-mbru-amr"
BEST_ROOT="$PROJECT/results/Best_Assemblies"
MANIFEST="$BEST_ROOT/master_manifest.tsv"

# Overwrite the old manifest
printf "Sample\tPlatform\tSpecies\tGroup\tFasta\tSource_QC_TSV\n" > "$MANIFEST"

add_group() {
    local platform="$1"
    local species="$2"

    local group="${platform}_${species}"
    local fasta_dir="$BEST_ROOT/$platform/$species"
    local source_tsv="$BEST_ROOT/${group}_selected.tsv"

    echo "Processing group: $group"
    echo "  FASTA directory: $fasta_dir"
    echo "  Source TSV:       $source_tsv"

    if [ ! -d "$fasta_dir" ]; then
        echo "WARNING: Directory not found, skipping: $fasta_dir" >&2
        return
    fi

    if [ ! -f "$source_tsv" ]; then
        echo "WARNING: Source QC TSV not found: $source_tsv" >&2
    fi

    find "$fasta_dir" \
        -maxdepth 1 \
        -type f \
        -name "*_contigs.fasta" \
        -print0 |
    sort -z |
    while IFS= read -r -d '' fasta
    do
        filename=$(basename "$fasta")
        sample="${filename%_contigs.fasta}"

        printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
            "$sample" \
            "$platform" \
            "$species" \
            "$group" \
            "$fasta" \
            "$source_tsv" \
            >> "$MANIFEST"
    done
}

add_group "ONT" "ecoli"
add_group "ONT" "kpneumoniae"
add_group "PacBio" "ecoli"
add_group "PacBio" "kpneumoniae"

echo
echo "Manifest created:"
echo "$MANIFEST"

echo
echo "Rows per group:"
awk -F'\t' '
NR > 1 {
    count[$4]++
}
END {
    for (group in count) {
        printf "%-25s %d\n", group, count[group]
    }
}
' "$MANIFEST" | sort

echo
echo "Total manifest rows:"
awk 'END {print NR-1}' "$MANIFEST"

echo
echo "Checking FASTA paths..."

missing=0

while IFS=$'\t' read -r sample platform species group fasta source_tsv
do
    if [ "$sample" = "Sample" ]; then
        continue
    fi

    if [ ! -f "$fasta" ]; then
        echo "MISSING FASTA: $fasta"
        missing=$((missing + 1))
    fi
done < "$MANIFEST"

if [ "$missing" -eq 0 ]; then
    echo "All FASTA paths exist."
else
    echo "ERROR: $missing FASTA file(s) are missing."
    exit 1
fi

echo
echo "First rows:"
head -n 6 "$MANIFEST"
