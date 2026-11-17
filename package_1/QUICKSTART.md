# Quickstart: install and run the continent and Archean examples

For anyone who just wants to build the code and reproduce the continent
formation (M1) and Archean verification runs.

## 1. Requirements

- Linux, gcc, and an MPI stack (e.g. OpenMPI: mpirun on the PATH)
- Python 3 with: numpy, scipy, matplotlib, mpi4py
  (pip install numpy scipy matplotlib mpi4py)

## 2. Build

From this directory (package_1/):

    gcc -O3 -march=native -ffast-math -fPIC -shared TLBM_fast.c -o TLBM_fast.so -lm

## 3. Run

Always from this directory, with an unlimited stack:

    ulimit -s unlimited

Archean single phase (512 x 64 cells, 600 My; hours on 16 cores):

    mpirun -np 16 python3 TLBM_fast.py TEST archean_cart_pre_stream

Archean two phase (same size; needs its initial-condition files first):

    mkdir -p projects/TEST/archean_cart_rk_shifted_redtop
    cp projects/TEST/rho_r_150km.npy projects/TEST/rho_b_150km.npy \
       projects/TEST/archean_cart_rk_shifted_redtop/
    mpirun -np 16 python3 TLBM_fast.py TEST archean_cart_rk_shifted_redtop

Continent formation M1 (5096 x 512 cell annulus, 3 Gyr; this is the large
run of the paper: allow several days, originally about 81 hours on an HPC
node; more ranks help):

    mpirun -np 16 python3 TLBM_fast.py PEPI_melt M1

## 4. Where the output goes

For a run <proj>/<run>:

- projects/<proj>/<run>.out      log (progress, mach, diagnostics)
- projects/<proj>/<run>/plots/   field images (png), written every Inc_plots
- projects/<proj>/<run>/data/    field snapshots (npy), every Inc_data

Early sanity check: the .out should show "Run is starting", then periodic
"Time = ..." lines with small mach values (about 0.007 at 100 My for M1).

Copyright (c) 2026 Peter Mora. See ../COPYRIGHT.txt / ../LICENSE.md.
