#!/bin/csh -f

set homedir = $homedir
set PROJECTS = examples
set RUNS     = "archean hadean mantle mantle_annulus mantle3D"
set VARS     = "T u"
set fps      = `cat $homedir/FPS`

if(($#argv == 1 && $1 == "-h") || $#argv > 1) then

    echo "------------------------------------------------------------------------"
    echo "Make movies from plot files of all examples runs"
    echo " "
    echo "    project = $PROJECTS"
    echo "    runs    = ${RUNS}"
    echo "    vars    = ${VARS}"
    echo "------------------------------------------------------------------------"
    echo "Usage: movies_ex [fps]"
    echo "------------------------------------------------------------------------"
    echo "Examples: movies_ex"
    echo "          movies_ex 20"
    echo "------------------------------------------------------------------------"
    echo "Notes: (1) fps defaults to the FPS file value of ${fps} frames per second if unspecified"
    echo "       (2) Movies will be created for variables in file VARS (eg. T, u, etc)"
    echo "           using the plot files in the plots directory, ie. T_0*.png etc"
    echo "       (3) Plot  files (*.png) must exist      : homedir/projects/project/run/plots/"
    echo "       (4) Movie files (*.avi) will be created : homedir/projects/project/run/movies/"
    echo "       (5) Runs can be done using the run alias or by running GTXLBM_run.py"
    echo "       (6) Plots files can be created from data files using the plots alias or by running plot_npy.py"
    echo "------------------------------------------------------------------------"
    exit
endif

if($#argv == 1) then
    set fps = $1
endif

echo "========================================================================"
echo "Generating movies  at ${fps} frames per second for projects specified in"
echo "file PROJECTS = [${PROJECTS}], and for all runs specified in file"
echo "RUNS = [${RUNS}]"
echo "for variables specified in file VARS = [${VARS}]"

set next_project = false
foreach project0 ($PROJECTS)
    set add_run      = false
    set in_out_added = false
    echo "========================================================================"
    echo "Generating movies for project = $project0"
    foreach run ($RUNS)
        if($run != $project0) then
            set project = $project0
            set next_project = true
            set add_run      = true
            set next_run     = false
        endif
        if($add_run == true && $run != $project0) then
            set run_dir = ./projects/$project/$run/
            if($run != $project0) then
                set run_dir = ./projects/$project/$run/
                if(-e $run_dir) then
                    if($project == examples & $run == mantle3D) then
                        movie3D.s xy 7,30,53 $fps $run $project
                    else
                        movie2D.s $fps $run $project
                    endif
                endif
            endif
            set next_run     = true
        endif
        set next_project = false
    end
end
echo "============================================================================"
echo "All movies have been made!"
echo "============================================================================"
