#!/bin/csh -f

set PROJECTS = `cat $homedir/PROJECTS`
set RUNS     = `cat $homedir/RUNS`

if($#argv > 0) then

    echo "-------------------------------------------------------------------"
    echo "Generate plot files from data files for multiple run(s) & project(s)"
    echo "-------------------------------------------------------------------"
    echo "Usage:   plots"
    echo "-------------------------------------------------------------------"
    echo "Notes: (1) Run(s) and project(s) per RUNS & PROJECTS files"
    echo "       (2) Data files (*.npy) must exist      in $homedir/projects/project/run/data/"
    echo "       (3) Plot files (*.png) will be created in $homedir/projects/project/run/plots/"
    echo "-------------------------------------------------------------------"
    exit
endif

echo "========================================================================"
echo "Submitting background jobs to make plots for projects specified in file"
echo "PROJECTS = [${PROJECTS}], and for all runs specified in file RUNS ="
echo "[${RUNS}]"

pushd $homedir >& /dev/null

set next_project = false
foreach project0 ($PROJECTS)
    set add_run      = false
    set in_out_added = false
    echo "========================================================================"
    echo "Submitting jobs to generate plots for project = $project0"
    echo "------------------------------------------------------------------------"
    foreach run ($RUNS)
        if($run != $project0) then
            set project = $project0
            set next_project = true
            set add_run      = true
        endif
        if($add_run == true && $run != $project0) then
            set run_dir = ./projects/$project/$run/
            if($run != $project0) then
                set run_dir = ./projects/$project/$run/
                if(-e $run_dir) then
                    echo "Submitting job to generate plots for run = $run"
                   python3 plot_npy.py $project $run >& ./projects/$project/$run/plots/plots.out &
                else
                    echo "Run directory does not exist for run = $run, ie. $run_dir does not exist!"
                endif
            endif
        endif
        set next_project = false
    end
end
echo "========================================================================"
echo "Background jobs to make plots for projects specified in file PROJECTS ="
echo "[${PROJECTS}], and for all runs specified in file RUNS ="
echo "[${RUNS}]"
echo "------------------------------------------------------------------------"
echo "For the status of any given plot generation run,"
echo "    type: tail -20 $homedir/projects/project/run/plots/plots.out"
echo      or:   use the statusp alias: statusp lines run project"
echo "========================================================================"

popd >& /dev/null
