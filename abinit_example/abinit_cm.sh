#!/bin/bash
#SBATCH --job-name="abinit-cm"
#SBATCH --time=01:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=1G
#SBATCH --partition=compute
#SBATCH --account=innovation

# Load modules:
module load 2025
module load miniconda3

# Activate conda, run job, deactivate conda
conda activate carriermult

date

srun python ../cmscript.py -a bulk_1o_GSR.nc

date
conda deactivate

