#!/bin/bash
#SBATCH --job-name=resample-stereo-fixed
#SBATCH --account=chemistry
#SBATCH --partition=kemi1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=01:00:00
#SBATCH --output=experiments/stereo-fixed-%j.out

set -euo pipefail
cd /lustre/hpc/kemi/hteahan/rxflow
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
exec .venv/bin/python -u experiments/check_rxembed_stereo.py \
    --output experiments/rxembed_730eae5_stereo_fixed --keep 5
