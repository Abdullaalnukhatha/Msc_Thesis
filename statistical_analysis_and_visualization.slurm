#!/bin/bash
#SBATCH --job-name=stat_analysis
#SBATCH --time=04:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --output=/rds/projects/e/elhamsak-mbru-amr/tempdata/statistical_analysis_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-mbru-amr/tempdata/statistical_analysis_%j.err

set -euo pipefail

PROJ="/rds/projects/e/elhamsak-mbru-amr"

ENV_DIR="$PROJ/envs/ref_valid"
PYTHON_BIN="$ENV_DIR/bin/python"
JUPYTER_BIN="$ENV_DIR/bin/jupyter"

NOTEBOOK="$PROJ/Accepted_scripts/Statistical_analysis_and_visualization.ipynb"
METRICS_FILE="$PROJ/results/plasmid_benchmark/validation_metrics/all_metrics.tsv"
OUTPUT_DIR="$PROJ/results/plasmid_benchmark/statistical_analysis"

mkdir -p "$PROJ/tempdata"

echo "============================================================"
echo " Statistical analysis and visualization"
echo "============================================================"
echo "Project      : $PROJ"
echo "Environment  : $ENV_DIR"
echo "Python       : $PYTHON_BIN"
echo "Jupyter      : $JUPYTER_BIN"
echo "Notebook     : $NOTEBOOK"
echo "Metrics file : $METRICS_FILE"
echo "Output dir   : $OUTPUT_DIR"
echo "Started      : $(date)"
echo "Host         : $(hostname)"
echo "============================================================"
echo ""

module purge
module load bluebear
module load bear-apps/2024a/live
module load Miniforge3/25.3.0-3

if [ ! -x "$PYTHON_BIN" ]; then
    echo "ERROR: Python not found: $PYTHON_BIN"
    exit 1
fi

if [ ! -x "$JUPYTER_BIN" ]; then
    echo "ERROR: Jupyter not found: $JUPYTER_BIN"
    exit 1
fi

if [ ! -f "$NOTEBOOK" ]; then
    echo "ERROR: Notebook not found: $NOTEBOOK"
    exit 1
fi

if [ ! -f "$METRICS_FILE" ]; then
    echo "ERROR: Metrics file not found: $METRICS_FILE"
    exit 1
fi

echo "Checking required Python packages..."
"$PYTHON_BIN" - <<'PY'
packages = [
    "pandas",
    "numpy",
    "matplotlib",
    "seaborn",
    "scipy",
    "statsmodels",
    "scikit_posthocs",
    "nbconvert",
    "ipykernel",
]
missing = []
for pkg in packages:
    try:
        m = __import__(pkg)
        print("OK:", pkg, getattr(m, "__version__", "version not shown"))
    except Exception as e:
        print("MISSING:", pkg, e)
        missing.append(pkg)

if missing:
    raise SystemExit("Missing required packages: " + ", ".join(missing))
PY

echo ""
echo "Replacing statistical analysis output folder:"
echo "$OUTPUT_DIR"
rm -rf "$OUTPUT_DIR"
mkdir -p "$OUTPUT_DIR"

EXECUTED_NOTEBOOK="$OUTPUT_DIR/Statistical_analysis_and_visualization_executed.ipynb"

echo ""
echo "Running notebook..."
"$JUPYTER_BIN" nbconvert \
    --to notebook \
    --execute "$NOTEBOOK" \
    --output "$EXECUTED_NOTEBOOK" \
    --ExecutePreprocessor.timeout=-1 \
    --ExecutePreprocessor.kernel_name=python3

echo ""
echo "============================================================"
echo " Statistical analysis complete"
echo "============================================================"
echo "Executed notebook : $EXECUTED_NOTEBOOK"
echo "Output folder     : $OUTPUT_DIR"
echo "Finished          : $(date)"
echo ""
echo "Generated files:"
find "$OUTPUT_DIR" -maxdepth 1 -type f -printf "  %f\n" | sort
