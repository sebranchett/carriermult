#!/bin/bash
#SBATCH --job-name="qe-cm"
#SBATCH --time=01:00:00
#SBATCH --ntasks=1  # this can only be 1 for cm script
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --partition=compute
#SBATCH --account=innovation

# Activate uv environment at head of repository
source ${PWD}/../.venv/bin/activate

date

srun python ../cmscript.py

date
