#!/bin/bash
#SBATCH --job-name=tlbm_ladder
#SBATCH --partition=cpu_x440
#SBATCH --nodes=1
#SBATCH --exclusive
#SBATCH --ntasks=64
#SBATCH --cpus-per-task=1
#SBATCH --time=04:00:00
#SBATCH --output=ladder_%j.log
source /etc/profile.d/modules.sh
module purge
module load conda/25.08 openmpi/gnu/5.0.10
cd $HOME/src/GTXLBM_BENCH
ulimit -s unlimited
gcc -O3 -ffast-math -fPIC -shared TLBM_fast.c -o TLBM_fast.so -lm
for NP in 4 8 16 32 64; do
  for CODE in fast prefuse; do
    R=bench_${CODE}_${NP}
    cp projects/PEPI_melt/M1_hpcb.in projects/PEPI_melt/$R.in
    echo "=== $R start $(date +%T) ==="
    mpirun -np $NP python3 TLBM_${CODE}.py PEPI_melt $R > /dev/null 2>&1
    grep -aE "Cores used|Mflups|Edges:" projects/PEPI_melt/$R.out 2>/dev/null | tail -6
    rm -rf projects/PEPI_melt/$R projects/PEPI_melt/$R.bmk projects/PEPI_melt/$R.growth
  done
done
echo "LADDER DONE"
