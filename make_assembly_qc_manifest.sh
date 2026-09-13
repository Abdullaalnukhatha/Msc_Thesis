#!/bin/bash

set -euo pipefail
export LC_ALL=C

# ============================================================
# Create the one-column manifest required by:
# assembly_qc_array.slurm
#
# Each manifest line contains the absolute path to one non-empty
# assembly FASTA file.
#
# Expected assembly layouts include:
#
# results/assemblies/ONT/ecoli/SRRxxxx_contigs.fasta
# results/assemblies/ONT/kpneumoniae/SRRxxxx_contigs.fasta
# results/assemblies/PacBio/ecoli/SRRxxxx_contigs.fasta
# results/assemblies/PacBio/kpneumoniae/SRRxxxx_contigs.fasta
#
# The reference FASTA is NOT placed in this manifest.
# The QC script looks it up separately using:
# data/ReferenceFASTA/reference_mapping.tsv
# ============================================================

PROJ="/rds/projects/e/elhamsak-mbru-amr"
ASSEMBLY_DIR="$PROJ/results/assemblies"
TEMP_DIR="$PROJ/tempdata"
MANIFEST="$TEMP_DIR/assembly_qc_manifest.tsv"
REFERENCE_MAP="$PROJ/data/ReferenceFASTA/reference_mapping.tsv"

mkdir -p "$TEMP_DIR"

# ------------------------- Validation -------------------------

[[ -d "$ASSEMBLY_DIR" ]] || {
    echo "ERROR: Assembly directory does not exist:"
    echo "  $ASSEMBLY_DIR"
    exit 1
}

[[ -s "$REFERENCE_MAP" ]] || {
    echo "ERROR: Reference mapping file is missing or empty:"
    echo "  $REFERENCE_MAP"
    exit 1
}

# Write to a temporary file first so an existing valid manifest is
# not destroyed if manifest creation fails.
TMP_MANIFEST="${MANIFEST}.tmp.$$"

trap 'rm -f "$TMP_MANIFEST"' EXIT

# ---------------------- Create manifest -----------------------

find "$ASSEMBLY_DIR" \
    -type f \
    -name "*_contigs.fasta" \
    -size +0c \
    -print \
    | sort -u \
    > "$TMP_MANIFEST"

TOTAL="$(wc -l < "$TMP_MANIFEST")"
TOTAL="${TOTAL//[[:space:]]/}"

if (( TOTAL == 0 )); then
    echo "ERROR: No non-empty *_contigs.fasta files were found under:"
    echo "  $ASSEMBLY_DIR"
    exit 1
fi

# Replace the previous manifest only after successful creation.
mv "$TMP_MANIFEST" "$MANIFEST"
trap - EXIT

# -------------------------- Summary ---------------------------

echo "============================================================"
echo "Assembly QC manifest created successfully"
echo "============================================================"
echo "Assembly directory : $ASSEMBLY_DIR"
echo "Manifest           : $MANIFEST"
echo "Total assemblies   : $TOTAL"
echo

echo "Assemblies by platform and species:"

awk -F '/' '
{
    platform = $(NF-2)
    species  = $(NF-1)
    count[platform "\t" species]++
}
END {
    for (group in count) {
        print group "\t" count[group]
    }
}' "$MANIFEST" \
    | sort \
    | awk -F '\t' 'BEGIN {
        printf "%-12s %-18s %s\n", "PLATFORM", "SPECIES", "COUNT"
        printf "%-12s %-18s %s\n", "--------", "-------", "-----"
    }
    {
        printf "%-12s %-18s %d\n", $1, $2, $3
    }'

echo
echo "First five manifest entries:"
head -n 5 "$MANIFEST"

echo
echo "Last five manifest entries:"
tail -n 5 "$MANIFEST"

echo
echo "Reference mapping file:"
echo "  $REFERENCE_MAP"
echo
echo "The manifest contains assembly paths only."
echo "The assembly QC array script will obtain the reference accession"
echo "from reference_mapping.tsv and locate the corresponding FASTA."
echo
echo "Submit the QC array with:"
echo
echo "  TOTAL=\$(wc -l < \"$MANIFEST\")"
echo "  sbatch --array=1-\${TOTAL}%20 assembly_qc_array.slurm"
echo "============================================================"
