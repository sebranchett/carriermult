#!/bin/bash

#SBATCH --job-name=abinit-dft
#SBATCH --partition=compute
#SBATCH --account=innovation
#SBATCH --time=01:00:00
#SBATCH --ntasks=4
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=3900MB

# find your account with:
# sacctmgr list -sp user $USER

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module use $(find /projects/electronic_structure/software/spack -name "Core")
module load 2025 python abinit

# Make sure you have the correct python packages installed for your
# python version. Only need to do this once:
# module load python
# python -m pip install numpy pyyaml pandas --user

export ABI_HOME=$ABINIT_ROOT
export ABI_TESTS=$ABI_HOME/tests/
export ABI_PSPDIR=${PWD}/NCFR_Pseudos/

date

time srun abinit < x.files > log 2> err

date

