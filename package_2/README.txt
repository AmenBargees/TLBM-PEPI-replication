============================================================================================
OVERVIEW OF GTXLBM FOR MODELLING FLUID DYNAMICS
============================================================================================

Requirements: Linux operating system with mpiexec plus Python and Python libraries: numpy, math,
              random, scipy, matplotlib, ctypes, mpi4py, os, time, sys and datetime. If mpiexec,
              python3 or Python libraries are not installed, install them via apt such as:

                  sudo apt update
                  sudo apt install -y openmpi-common libopenmpi-dev (install mpi)
                  sudo apt install python3                          (install python3)
                  sudo apt install python3-mpi4py                   (install mpi4py library)
                  ...

               GTXLBM uses scripts written in the c-shell to run, display plots and more
               so if csh is not installed, install it via: sudo apt install csh

Installation: (1) Extract GTXLBM files from the tar file such as: tar xvf GTXLBM_V2.96_Peter_Mora_2026-08-17.tar

              (2) Type the command "install.s". This will create GTXLBM_run.py - the Python code used
                  for runs - and compile the GTXLBM C kernel using the gcc compiler - the GNU C
                  compiler - which is standard on Linix systems. If your computer is equipped with
                  Intel chips, optimal performance requires the codename of the Intel chip to specified
                  to the gcc compiler. Typing the command "make.s -h" will display the help documentation
                  for make.s which explains how to determine what the chip on your computer and its
                  codename. Once the codename is known, the command "make.s codename" to compile the
                  GTXLBM C kernel and achieve the best performance on your computer.

               (3) GTXLBM_run.py is now ready to be run, but first, it is important to understand
                   the drectory

Running:      There are three ways to do GTLBMX runs as below but these all require a run input file and a project to exist.
              A run input file named run_name.in must exist in a project directory named ./projects/project_name/ in order to do a run.
              Examples of run inputt files can be found in ./projects/examples.

              (1) Via running GTXLBM_run.py using mpiexec

                  mpiexec -np cores python3 GTXLBM_run.py project_name run_name

                  where:      cores        = number of cores
                              project_name = project name
                              run_name     = run     name

              (2) Via the run script where the run and project names are per the RUN & PROJECT files, and the number of cores is per the CORES file

                  run

              (3) Via the runs script where the run and project name(s) are per the RUNS & PROJECTS files, and the number of cores is per the CORES file

                  runs

Output:        GTXLBM runs generate an output file showing the input parameters that were read from the input file,
              initializations, and regular updates on run status such as the physical time modelled so far, the run time so far
              and total run time to completion, the mach number (which relates to accuracy) far, and the output plot (.png) and
              data (.npy) files. Plot and data files at specified time increments (eg. 10 My) are written to directories
              named ./projects/project_name/run_name/plots and ./projects/project_name/run_name/data, and have names
              such as T_0000000500_My.png and u_0000000500_My.png for temperature and velocity at 500_My.

Plots:        Saving plot files can slow GTXLBM so it is best to save only enough plot files to check run progress and/or
              to create a basic movie. Saving data files is fast and can be done more frequently. As such, it is best to save
              ony of order 30 - 50 plot files during the run and to generate more plot files from the saved data files should
              this be required. This can be achieved by running plot_npy.py in three ways similarly to GTXLBM_run.py as below.

              (1) python3 plot_npy.py project_name run_name                    (for 2D runs)
                  python3 plot_npy.py project_name run_name coordinates        (for 3D runs where the plane is assumed to be xy)
                  python3 plot_npy.py project_name run_name plane coordinates  (for 3D runs where the plane can be xy, xz or yz)

              (2) Via the plot script where the run and project names are per the RUN & PROJECT files, and the number of cores is per the CORES file

                  plot

              (3) Via the plots script where the run and project name(s) are per the RUNS & PROJECTS files, and the number of cores is per the CORES file

                  plots

Scripts:      There many other useful scripts named scriptname.s located in directory ./scripts, each with an alias of script_name.
              These scripts and aliases make it easy to do runs, check run status, display plots, generate movie files from plot
              files (ie. generate .avi files in the movie directory ./projects/project_name/run_name/movies/), pushd to the plots,
              data, movies or project directories and more! The scripts explained above were run, runs, plot and plots, and others include:

              movie  (make movie files from plot files for variables specified in the VARS file for the run & project specified in
                        files RUN and PROJECT at the number of frames per second specified in the FPS file)
              movies (same as movie but for the multiple runs. Namely, the run(s) & project(s) specified in files RUNS and PROJECTS)

              ep plot_file.png (display plot_file.png in the plots directory for the run & project specified in files RUN and PROJECT via eog)
                               (useful as eog does not require an X-server)
              dp plot_file.png (display plot_file.png in the plots directory for the run & project specified in files RUN and PROJECT via display)
              lsp              (use ls to list the plots directory for the run & project specified in files RUN and PROJECT)
              sr               (tail the output file of the run & project specified in files RUN and PROJECT to get the run status)
              sp               (tail the output file of a plot generation run using plot_npy.py for the run & project specified in files RUN and PROJECT to the status)
              prr              (more the output file for the run & project specified in files RUN and PROJECT)

              tar.s    (generate a tarfile such as GTXLBM_V2.96_Peter_Mora_2026-08-17_username_YYYY-MM-DD.tar in directory ./tarfiles/ where
                       where 2.95 is the GTXLBM version, Peter_Mora is the GTXLBM author, username is your username on your computer, and
                       YYYY-MM-DD specify the date you created this tarfile. The tarfile will contain the source directory where you installed
                       GTXLBM, along with the input and output files for the run(s) and project(s) specified in files RUNS and PROJECTS, plus
                       the plots, data or movies directories specified in file DIRS for these runs(s) & projects(s). Useful to save
                       run output, plots and data to collaborators, or for long termm backup of projects.

              bak      (generate a backup directory such as ./backups/GTXLBM_V2.95_Peter_Mora_2026-08-17_username_YYYY-MM-DD/
                        This is the same as tar.s but the tarfile is created and then extracted in the backup directory which allows
                        instant access to files in the backup directory.
                        
Documentation: All GTXLBM programs and scripts are self-documenting when passed -h as a sole argument. For example,

               mpiexec -np 1 python3 GTXLBM_run.py -h
               python3 plot_npy.py -h
               plot.s  -h   (note: the -h help argument will work for all scripts (scriptname.s), but not always via their
                             aliases as some aliases source the script file. For example, plot -h will fail!

Welcome and good luck!

Peter Mora

============================================================================================
Example runs and input files
============================================================================================

The examples below can all be run using the run_examples script.

archean.in:        Archean geodynanmics with temperature dependent viscosity and a strongly nonlinear rheology showing
                   the transition to plate tectonics: stagnant lid convection - crustal thinning - drip tectonics - plate tectonics.

hadean.in          Hadean geodynamics to illustrate modelling high Rayleigh numbers in the turbulent regime.

mantle.in          2D solid-state mantle convection in a full annulus.
mantle_annulus.in  2D solid-state mantle convection in a part annulus spanning 120 degrees.
mantle3D.in        3D whole mantle convection in a Cartesian model.

============================================================================================
DEVELOPERS
============================================================================================

See README_DEVELOPERS.txt

============================================================================================
UNDERSTANDING THE TLBM:
============================================================================================

See papers in directory ./papers/. The file PAPERS_LIST.txt contains a concise
overview of what each paper covers.

============================================================================================
