#!/bin/csh -f

if($#argv == 1 && $1 == "-h" || $#argv > 3) then

    echo "-------------------------------------------------------------------"
    echo "Background run of GTXLBM"
    echo "-------------------------------------------------------------------"
    echo "Usage: run [cores run project]"
    echo "-------------------------------------------------------------------"
    echo "Examples: run"
    echo "          run archean"
    echo "          run mantle 32"
    echo "          run mantle3D 80 examples"
    echo "-------------------------------------------------------------------"
    echo "Notes: (1) Run & project per RUN and PROJECT files if unspecified"
    echo "       (2) Cores per CORES file if unspecified"
    echo "       (3) Type tr or tr.s to tail the .out file & check on run status"
    echo         (4) Run & project per RUN & PROJECT files if unspecified"
    echo "       (5) Requires the run input file: ./projects/project/run.in"
    echo "       (5) Output files & directories :"
    echo " "
    echo "           ./projects/project/run.out"
    echo "           ./projects/project/run/data/*.npy"
    echo "           ./projects/project/run/plots/*.png"
    echo "           ./projects/project/run/movies/"
    echo "------------------------------------------------------------------------"
    exit
endif

set homedir = $homedir
set cores   = `cat $homedir/CORES`
set run     = `cat $homedir/RUN`
set project = `cat $homedir/PROJECT`

if ($#argv >= 1) then
    set cores = $1
endif
if($#argv >= 2) then
    set run   = $2
endif
if($#argv == 3) then
    set project = $3
endif

pushd $homedir >& /dev/null

set run_in_file = $homedir/projects/$project/$run.in
set project_dir = $homedir/projects/$project/

echo "------------------------------------------------------------------------"
if (!(-e $project_dir)) then
    echo "ERROR: Project directory does not exist!"
    echo "           $project_dir"
    echo "------------------------------------------------------------------------"
endif
if(!(-e $run_in_file)) then
    echo "ERROR: Input file does not exist!"
    echo "       Run could not be done for non-existent input file:"
    echo "           $run_in_file"
else
    echo "Input file exists: ${run_in_file}"
    echo "Preparing to start run = $run, project = $project ..."
    echo "    File & directory cleanup ..."
    echo -n "        rm $homedir/projects/$project/$run.out ... "
    rm $homedir/projects/$project/$run.out        >& /dev/null
    echo "Done"

    echo -n "        rm $homedir/projects/$project/$run.bmk ... "
    rm $homedir/projects/$project/$run.bmk        >& /dev/null
    echo "Done"

    echo -n "        rm -r $homedir/projects/$project/$run/plots/ ... "
    rm -r $homedir/projects/$project/$run/plots/  >& /dev/null
    echo "Done"

    echo -n "        rm -r $homedir/projects/$project/$run/data/ ... "
    rm -r $homedir/projects/$project/$run/data/   >& /dev/null
    echo "Done"

    echo -n "        rm -r $homedir/projects/$project/$run/movies/ ... "
    rm -r $homedir/projects/$project/$run/movies/ >& /dev/null
    echo "Done"

    echo -n "Starting run on $cores cores ..."
    date
    echo "Output file: $homedir/projects/$project/$run.out"
    echo "mpiexec -np $cores python3 GTXLBM_run.py $project $run >& $homedir/projects/$project/$run.err & ... running ..."
    mpiexec -np $cores python3 GTXLBM_run.py $project $run >& $homedir/projects/$project/$run.err &
endif
echo "-------------------------------------------------------------------"

popd >& /dev/null
