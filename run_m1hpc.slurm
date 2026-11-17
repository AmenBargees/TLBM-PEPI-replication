#!/bin/bash
#SBATCH --job-name=M1_hpc
#SBATCH --partition=cpu_x440
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --output=M1_hpc_%j.log
source /etc/profile.d/modules.sh
module purge
module load conda/25.08 openmpi/gnu/5.0.10
cd $HOME/src/GTXLBM_BENCH
ulimit -s unlimited
gcc -O3 -ffast-math -fPIC -shared TLBM_fast.c -o TLBM_fast.so -lm
mpirun -np 64 python3 TLBM_fast.py PEPI_melt M1_hpc
