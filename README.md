# PEPI paper replication package

One code tree (TLBM V2.6), compiled with:
  gcc -O3 -march=native -ffast-math -fPIC -shared TLBM_fast.c -o TLBM_fast.so -lm
Run pattern (16 MPI ranks, as in the paper):
  ulimit -s unlimited && mpirun -np 16 python3 TLBM_fast.py <project> <run>

## TLBM_V26/
- projects/RT/vk_b5L.in                           van Keken immiscible experiment
- projects/RT/van-keken-miscible-2x.in            van Keken miscible experiment
- projects/PEPI_melt/M1.in                        Archean continent formation example
- projects/TEST/archean_cart_pre_stream.in        single phase archean
- projects/TEST/archean_cart_rk_shifted_redtop.in two phase archean


## Verification record (2026-11, KFUPM HPC, 16 ranks)

Each input was rerun under this single code version and checked against the
original run's output:

- vk_b5L:                   Vrms 0.000753 at 500 Ky (original 0.00075)
- van-keken-miscible-2x:    mach 0.00499 at 25 Ky (original 0.005)
- M1:                       mach 0.00672 / 0.0134 at 100 / 200 Ky
                            (original 0.0067 / 0.013)
- archean_cart_rk_shifted_redtop: rho_r field at 90 My matches the paper
                            figure blob for blob
- archean_cart_pre_stream:  statistically equivalent (random boundary
                            initialization draws a different realization;
                            cluster and PC31 agree bit for bit with each other)

archean_cart_rk_shifted_redtop requires rho_r_150km.npy and rho_b_150km.npy
(included in projects/TEST/) copied into its run directory
projects/TEST/archean_cart_rk_shifted_redtop/ before launch.

Copyright (c) 2026 Peter Mora. See COPYRIGHT.txt / LICENSE.md.
