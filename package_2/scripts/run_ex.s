#!/bin/csh -f

set PROJECT  = examples
set RUNS     = "archean hadean mantle mantle_annulus mantle3D"

if(($#argv == 1 && $1 == "-h") || $#argv > 1) then
    echo "------------------------------------------------------------------------"
    echo "Run GTXLBM example runs  : $RUNS"
    echo "           output project: $PROJECT" 
    echo "------------------------------------------------------------------------"
    echo "Usage:    run_ex [cores]"
    echo "------------------------------------------------------------------------"
    echo "Examples: run_ex"
    echo "          run_ex 16"
    echo "------------------------------------------------------------------------"
    echo "Notes: (1) The number of cores used is up to 80 if the number of cores"
    echo "       (2) GTXLBM will create output files & directories:"
    echo "               ./projects/project/run.out"
    echo "               ./projects/project/run/data/*.npy"
    echo "               ./projects/project/run/plots/*.png"
    echo "               ./projects/project/run/movies/"
    echo "       (3) run_ex is an alias to run_ex.s"
    echo "       (4) This script runs in the foreground sequentially so it is best"
    echo "           to run it in a separate terminal from the one you are working"
    echo "           on so you can continue to work on your current terminal."
    echo "       (5) Use the sr alias to check the status of a run"
    echo "------------------------------------------------------------------------"
    exit
endif
set RUNS     = "examples ${RUNS}"
echo "------------------------------------------------------------------------"
echo "About to start example runs using up to 80 cores by default. If your"
echo "computer has less than 80 cores, you must specify the number of cores"
echo "------------------------------------------------------------------------"
echo "Your workstation has cores as per the specs from lscpu below:"
lscpu | grep "Core(s) per socket: "
lscpu | grep "Socket(s): "
echo "------------------------------------------------------------------------"
echo "where cores = cores/socket x sockets"
echo "------------------------------------------------------------------------"
echo -n "Do you have enough cores to do the runs? "
set response = $<
if($response == n) then
    echo "------------------------------------------------------------------------------"
    echo ""
    echo "In that case, you must specify the number of cores up to the number you have"
    echo ""
    echo "Exiting!"
    echo "------------------------------------------------------------------------------"
    exit
endif

foreach RUN ($RUNS)
    if($RUN != examples) then
        if($#argv == 0) then
            if($RUN == archean) then
                set CORES = 16
            else if($RUN == hadean || $RUN == mantle || $RUN == mantle_annulus) then
                set CORES = 32
            else if($RUN == mantle3D) then
                set CORES = 80
            endif
        else if($#argv == 1) then
            set CORES = $1
        endif
        echo "============================================================"
        echo "Starting run $RUN of project $PROJECT on ${CORES} cores."
        echo "Use the sr alias see the run status or use the pr alias to"
        echo "more the output file $homedir/projects/$PROJECT/$RUN.out"
        echo "============================================================"
        cleanr.s $RUN $PROJECT
        echo 'Running "mpiexec -np '"${CORES} python3 GTXLBM.py $PROJECT $RUN"'"'
        mpiexec -np $CORES python3 GTXLBM.py $PROJECT $RUN
   endif
end
echo "============================================================"
echo "All example runs are complete!"
echo "============================================================"
