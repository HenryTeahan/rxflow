#!/bin/bash
#SBATCH --job-name=rxts-gxtb-92
#SBATCH --account=chemistry
#SBATCH --partition=kemi1
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --output=experiments/job_92_gxtb_test/slurm-%j.out

set -euo pipefail
cd /lustre/hpc/kemi/hteahan/rxflow
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec .venv/bin/python -u experiments/job_92_gxtb_test/run_test.py
