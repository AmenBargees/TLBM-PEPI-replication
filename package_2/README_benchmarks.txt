============================================================================================
CONCISE INSTRUCTIONS TO REPLICATE THE ANNULAR BENCHMARKS THE TLBM AGAINST ASPECT
============================================================================================

1. Type "install.s" on your linux computer to install GTXLBM.
2. Follow the instructions in README.txt to install any missing phython libraries and csh.
3. If you have less than 16 cores, edit the file CORES to the number of cores on your computer.
   If you have more than 16 cores, you may also edit the CORES file, but keep in mind that
   as many runs have small meshes, these runs will likely not run faster on more cores.
   On my 16 core workstation, all free-slip and zero-slip runs for Ra=10^4-10^9 take less
   than 12 hours.
4. From a c-shell, type: runs
5. After all of the benchmark runs have completed, type: python3 plot_bmk.py
   This generates the bmk plot files showing the scaling plots.
6. The benchmark scaling plots will be located in ./projects/benchmarks_annulus/bmk
7. The file README.txt contains more detailed instructions including an explanation
   of the file structure so you can locate plots of the input and output files of runs,
   and the temperature and velocity field png files produced in the runs.
8. The project and run names that will run are located in the file RUNS and PROJECTS, but
   if you change this file, see RUNS_benchmarks and PROJECTS_benchmarks.

============================================================================================
