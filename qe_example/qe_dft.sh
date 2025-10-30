#!/bin/bash

#SBATCH --job-name=qe-dft
#SBATCH --partition=compute
#SBATCH --account=innovation
#SBATCH --time=01:00:00
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3900MB
# #SBATCH --mail-type=ALL

# find your account with:
# sacctmgr list -sp user $USER

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

module use --append /projects/electronic_structure/.modulefiles
module load qe
module load yambo

WORKDIR=${PWD}/MoTe2
cd "$WORKDIR"

# Find the prefix name from scf.in
prefix=$(grep prefix scf.in | awk '{print $3}' | tr -d '[:punct:]')

# DFT with Quantum Espresso
# parameters to converge in combination with GW/BSE: ecutwfc, nbnd, K_POINTS
mkdir -p output
# scf
srun pw.x < scf.in > output/scf.out
# #SBATCH --mail-type=ALL

# Convert Quantum Espresso output to Yambo input
cd ${prefix}.save
srun p2y > ../output/p2y.out
cd ..
# p2y.out is empty if all went well
cp -rf ${prefix}.save/SAVE SAVE

# Create initialisation file (init.in) for Yambo and run the initialisation
rm -f init.in
yambo -i -V RL -F init.in  # this creates the input file
srun yambo -F init.in -J output/init  # this runs the initialisation
# Yambo report at output/r-init_setup

