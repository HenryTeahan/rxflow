#!/bin/bash
#SBATCH --job-name=rxflow
#SBATCH --array=0-99%100
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --output=logs/orca-%A_%a.out
#SBATCH --error=logs/orca-%A_%a.err

set -euo pipefail
# Submit from the repository root: sbatch scripts/hpc_worker.sbatch [config.toml]
# Ensure logs/ exists before submitting; Slurm opens the logs before this script runs.
# Load your cluster's ORCA/MPI modules and export ORCA_EXE,
# ORCA_EXTERNAL_TOOLS and XTB_EXE before submitting.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec .venv/bin/rxflow --config "${1:-hpc.toml}" worker
