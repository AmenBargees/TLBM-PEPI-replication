#!/bin/csh -f

if($#argv > 0) then

    echo "-------------------------------------------------------------------"
    echo "Run plot_bmk"
    echo "-------------------------------------------------------------------"
    echo "Usage: runb"
    echo "-------------------------------------------------------------------"
    echo "Example: runb"
    echo "-------------------------------------------------------------------"
    echo "Notes: (1) Runs & projects per specifications in plot_bmk.py"
    echo "       (2) runb is an alias to runb.s"
    echo "       (3) Output files & directories :"
    echo " "
    echo "           ./projects/project/bmk/*.png"
    echo "------------------------------------------------------------------------"
    exit
endif

set homedir  = $homedir
set runs     = `cat $homedir/RUNS`
set projects = `cat $homedir/PROJECTS`

pushd $homedir >& /dev/null

python3 plot_bmk.py

popd >& /dev/null
