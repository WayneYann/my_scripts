#!/bin/bash
#BSUB -J IHsolv[1-67]%10
#BSUB -o output.log
#BSUB -R "span[hosts=1]"
#BSUB -n 10
#BSUB -q west
#BSUB -e error_log

export GAUSS_SCRDIR="/scratch/slakman.b"
python solvation_energies.py
