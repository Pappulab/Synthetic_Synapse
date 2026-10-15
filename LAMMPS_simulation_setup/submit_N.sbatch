#!/bin/bash

#SBATCH -A engr-lab-pappu
#SBATCH -o LAMMPS_OUTPUT.%J
#SBATCH -J int0.1_het15.0
#SBATCH -n 1
#SBATCH -N 1
#SBATCH -p condo-pappu,general-cpu
#SBATCH --mem=8GB
#SBATCH --constraint="intel"

PROJECT_PATH=$PWD
echo "${PROJECT_PATH}"
cd "${PROJECT_PATH}"

module add engr-lab-pappu
module add intel-compiler-2020
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/project/fava/lib

CONFIG_FILE=lammps_sim.in
LAMMPS_PATH=/engrfs/project/pappu/gmitra/tests/lammps_6_29/build/lmp

# Run the program
"${LAMMPS_PATH}" -i "${CONFIG_FILE}"

