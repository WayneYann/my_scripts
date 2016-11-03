#!/bin/bash
#SBATCH --job-name=si-all
#SBATCH --array=1-90%5
#SBATCH -N 1
#SBATCH -n 20
#SBATCH --partition=west
#SBATCH --error=error_log_%a
#SBATCH --output=error_log_%a

export GAUSS_SCRDIR="/gss_gpfs_scratch/slakman.b"
python input.py
