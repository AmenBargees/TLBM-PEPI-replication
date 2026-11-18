#!/bin/csh -f

set PROJECTS = examples
set RUNS     = "examples archean hadean mantle mantle_annulus mantle3D"

if($#argv > 0) then
    echo "-------------------------------------------------------------------"
    echo "Generate plot files from data files of the example runs"
    echo "-------------------------------------------------------------------"
    echo "Usage:   plot_ex"
    echo "-------------------------------------------------------------------"
    echo "Example: plot_ex"
    echo "-------------------------------------------------------------------"
    echo "Note: Data files of the example runs of the example runs are generated"
    echo "by the run_examples alias to the run_examples.s script which runs"
    echo "all of the example runs"
    echo "-------------------------------------------------------------------"
    exit
endif

echo "========================================================================"
echo "Submitting background jobs to make plots for the examples project"

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
                    echo '"python3 plot_npy.py"'" $project $run >& ./projects/$project/$run/plots/plots.out &"
#                    python3 plot_npy.py $project $run >& ./projects/$project/$run/plots/plots.out &
                endif
            endif
        endif
        set next_project = false
    end
end
echo "========================================================================"
echo "Background jobs to make plots were submitted for the examples project"
echo "------------------------------------------------------------------------"
echo "See $homedir/projects/examples/run/plots/plots.out for the plot generation"
echo "status or type sp lines run examples, where lines is how many lines to tail"
echo "========================================================================"

popd >& /dev/null
