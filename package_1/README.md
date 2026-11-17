# A multiphase Thermal Lattice Boltzmann Method for modelling planetary dynamics PEPI paper replication package 

One code tree (TLBM V2.6), compiled with:
  gcc -O3 -march=native -ffast-math -fPIC -shared TLBM_fast.c -o TLBM_fast.so -lm
New users: see QUICKSTART.md for a minimal install-and-run guide.

Run pattern (16 MPI ranks, as in the paper):
  ulimit -s unlimited && mpirun -np 16 python3 TLBM_fast.py <project> <run>

## TLBM_V26/
- projects/RT/vk_b5L.in                           van Keken immiscible experiment
- projects/RT/van-keken-miscible-2x.in            van Keken miscible experiment
- projects/PEPI_melt/M1.in                        Archean continent formation example
- projects/TEST/archean_cart_pre_stream.in        single phase archean
- projects/TEST/archean_cart_rk_shifted_redtop.in two phase archean
- 
Copyright (c) 2026 Peter Mora. See COPYRIGHT.txt / LICENSE.md.
