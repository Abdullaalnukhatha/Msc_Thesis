#!/bin/bash

set -euo pipefail

# --- User Inputs ------------------------------------------------------------------------
read -rp "Input file (e.g. gcf_list.txt): " INPUTFILE
read -rp "Output directory [default: data/FASTA]: " OUTDIR
read -rp "Output ZIP filename [default: genomes.zip]: " ZIPNAME
read -rp "Max retries [default: 10]: " MAX_RETRIES
read -rp "Retry delay in seconds [default: 30]: " RETRY_DELAY

# Apply defaults if empty
OUTDIR=${OUTDIR:-data/FASTA}
ZIPNAME=${ZIPNAME:-genomes.zip}
MAX_RETRIES=${MAX_RETRIES:-10}
RETRY_DELAY=${RETRY_DELAY:-30}

ZIPFILE="${OUTDIR}/${ZIPNAME}"

# Validate input file
if [ ! -f "$INPUTFILE" ]; then
    echo "ERROR: Input file '$INPUTFILE' not found."
    exit 1
fi

mkdir -p "$OUTDIR"

# --- Helper functions ---------------------------------------------------------------------
run_with_retry() {
    local label="$1"; shift
    local attempt=1

    while [ $attempt -le $MAX_RETRIES ]; do
        echo "[$(date '+%H:%M:%S')] $label — attempt $attempt/$MAX_RETRIES..."
        if "$@"; then
            return 0
        fi

        echo "Failed. Waiting ${RETRY_DELAY}s before retry..."
        attempt=$((attempt + 1))
        sleep "$RETRY_DELAY"
    done

    echo "ERROR: $label failed after $MAX_RETRIES attempts."
    exit 1
}

# --- Step 1: Download dehydrated package --------------------------------------------------
echo ""
echo "[$(date '+%H:%M:%S')] Step 1/3 — Downloading dehydrated package..."

run_with_retry "Dehydrated download" \
    datasets download genome accession \
        --inputfile "$INPUTFILE" \
        --include genome \
        --dehydrated \
        --filename "$ZIPFILE"

# --- Step 2: Unzip -------------------------------------------------------------------------
echo ""
echo "[$(date '+%H:%M:%S')] Step 2/3 — Unzipping dehydrated package..."
unzip -q -o "$ZIPFILE" -d "$OUTDIR"

# --- Step 3: Rehydrate ---------------------------------------------------------------------
echo ""
echo "[$(date '+%H:%M:%S')] Step 3/3 — Rehydrating genome files (resumes if interrupted)..."

run_with_retry "Rehydration" \
    datasets rehydrate --directory "$OUTDIR"

echo ""
echo "[$(date '+%H:%M:%S')] Done! Genomes saved to: $OUTDIR"
