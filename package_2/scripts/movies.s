#!/bin/csh -f

set PROJECTS = `cat $homedir/PROJECTS`
set RUNS     = `cat $homedir/RUNS`
set VARS     = `cat $homedir/VARS`
set fps      = `cat $homedir/FPS`

if(($#argv == 1 && $1 == "-h") || $#argv > 1) then

    echo "------------------------------------------------------------------------"
    echo "Make movies from plot files for run(s) and project(s)"
    echo "------------------------------------------------------------------------"
    echo "Usage: movies [fps]"
    echo "------------------------------------------------------------------------"
    echo "Examples: movies"
    echo "          movies 20"
    echo "------------------------------------------------------------------------"
    echo "Notes: (1) Run(s) and project(s) per RUNS & PROJECTS files"
    echo "       (2) Variables per VARS file"
    echo "       (3) fps per FPS file if unspecified"
    echo "       (4) Plot  files (*.png) must exist      in $homedir/projects/project/run/plots/"
    echo "       (5) Movie files (*.avi) will be created in $homedir/projects/project/run/movies/"
    echo "       (6) Run the plots alias or run plot_npy.py to create plot files from data files"
    echo "------------------------------------------------------------------------"
    exit
endif

if($#argv == 1) then
    set fps = $1
endif

echo "========================================================================"
echo "Making movies at ${fps} frames/second for run(s), project(s) and var(s)"
echo "------------------------------------------------------------------------"
echo "Variables  : $VARS"
set i_run = 0
foreach i ($RUNS)
    if($i_run == 1) then
        echo "run(s)     : $i"
    else if($i_run > 1) then
        echo "             $i"
    endif
    @ i_run = $i_run + 1
end
set i_project = 0
foreach i ($PROJECTS)
    if($i_project == 0) then
        echo "project(s) : $i"
    else
        echo "             $i"
    endif
    @ i_project = $i_project + 1
end
echo "========================================================================"

set next_project = false
foreach project0 ($PROJECTS)
    set add_run      = false
    set in_out_added = false
    foreach run ($RUNS)
        if($run != $project0) then
            set project = $project0
            set next_project = true
            set add_run      = true
            set next_run     = false
        endif
        if($add_run == true && $run != $project0) then
            set run_dir = $homedir/projects/$project/$run/
            if($run != $project0) then
                set run_dir = $homedir/projects/$project/$run/
                if(-e $run_dir) then
                    if($project == examples & $run == mantle3D) then
                        movie3D.s xy 7,30,53  $fps $run $project
                    else
                        movie2D.s $fps $run $project true
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
