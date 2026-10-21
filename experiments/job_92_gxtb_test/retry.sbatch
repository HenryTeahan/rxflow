#!/bin/bash
#SBATCH --job-name=gxtb-peak-retry-92
#SBATCH --account=chemistry
#SBATCH --partition=kemi1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=24:00:00
#SBATCH --output=experiments/job_92_gxtb_test/retry-%j.out

set -euo pipefail
cd /lustre/hpc/kemi/hteahan/rxflow
source ./exports
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec .venv/bin/python -u experiments/job_92_gxtb_test/retry_job92.py
