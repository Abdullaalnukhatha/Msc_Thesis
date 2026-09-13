#!/bin/bash
#SBATCH --job-name=select_best_assemblies
#SBATCH --time=02:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --output=/rds/projects/e/elhamsak-mbru-amr/tempdata/select_best_assemblies_%j.out
#SBATCH --error=/rds/projects/e/elhamsak-mbru-amr/tempdata/select_best_assemblies_%j.err

module purge
module load bluebear
module load bear-apps/2024a/live

cd /rds/projects/e/elhamsak-mbru-amr

python Accepted_scripts/select_best_assemblies.py
