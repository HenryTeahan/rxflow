#!/bin/bash
#SBATCH --job-name=rxflow-sample
#SBATCH --array=0-12%13
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --output=sampling-%A_%a.out

set -euo pipefail
# Load the cluster's modules and export XTB_EXE before submitting.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec .venv/bin/rxflow --config "${1:-hpc.toml}" sample --ligand-index "$SLURM_ARRAY_TASK_ID"
