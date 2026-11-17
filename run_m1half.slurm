#!/bin/bash
#SBATCH --job-name=M1_half
#SBATCH --partition=cpu_x440
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00
#SBATCH --output=M1_half_%j.log
source /etc/profile.d/modules.sh
module purge
module load conda/25.08 openmpi/gnu/5.0.10
cd $HOME/src/GTXLBM_BENCH
ulimit -s unlimited
mpirun -np 64 python3 TLBM_fast.py PEPI_melt M1_half
