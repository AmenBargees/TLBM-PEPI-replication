# A multiphase Thermal Lattice Boltzmann Method for modelling planetary dynamics

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.24333821.svg)](https://doi.org/10.5281/zenodo.24333821)

Replication packages for the PEPI paper (Bargees et al.).

- package_1/ : TLBM V2.6 code and the five inputs replicating the paper's
               runs (van Keken immiscible and miscible, Archean continent
               formation M1, and the two archean verification runs).
               See package_1/README.md for the build line and run pattern.
- package_2/ : GTXLBM V2.96 (2026-11-18) code and the twelve annulus
               benchmark inputs (Ra 1e4-1e9, free-slip and zero-slip,
               Pr = 1000) for the TLBM vs ASPECT benchmark paper.
               Start with package_2/README_benchmarks.txt. Data and plot
               outputs are not included, only code and .in files.

Copyright (c) 2026 Peter Mora. See COPYRIGHT.txt / LICENSE.md.
