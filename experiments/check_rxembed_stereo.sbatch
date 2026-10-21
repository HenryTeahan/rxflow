#!/bin/bash
#SBATCH --job-name=rxembed-stereo-check
#SBATCH --account=chemistry
#SBATCH --partition=kemi1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:45:00
#SBATCH --output=experiments/rxembed-stereo-%j.out

set -euo pipefail
cd /lustre/hpc/kemi/hteahan/rxflow
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec .venv/bin/python -u experiments/check_rxembed_stereo.py
