#!/bin/bash
#SBATCH --job-name=rxflow
#SBATCH --array=0-99%100
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --output=hpc_results/slurm-%A_%a.out

set -euo pipefail
# Load your cluster's ORCA/MPI modules and export ORCA_EXE,
# OET_GXTB_ALPB_EXE and XTB_EXE before submitting.
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec .venv/bin/python hpc_workflow.py --config hpc.toml worker
