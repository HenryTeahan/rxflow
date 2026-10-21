#!/bin/bash
#SBATCH --job-name=rxflow-sample
#SBATCH --array=0-12%13
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --output=logs/sampling-%A_%a.out
#SBATCH --error=logs/sampling-%A_%a.err

set -euo pipefail
# Submit from the repository root: sbatch scripts/hpc_sample.sbatch [config.toml]
# Ensure logs/ exists before submitting; Slurm opens the logs before this script runs.
# Load the cluster's modules and export XTB_EXE before submitting.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec .venv/bin/rxflow --config "${1:-hpc.toml}" sample --ligand-index "$SLURM_ARRAY_TASK_ID"
