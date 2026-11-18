#!/bin/csh -f

set fps = `cat $homedir/FPS`

if(($#argv == 1 && $1 == "-h") || $#argv > 4) then
    echo "----------------------------------------------------------------------------"
    echo "Make movies from plot files of a run of a project of a 2D GTXLBM run"
    echo "----------------------------------------------------------------------------"
    echo "Usage:    movie.s [fps run project]"
    echo "----------------------------------------------------------------------------"
    echo "Examples: movie"
    echo "          movie 5"
    echo "          movie 10 archean"
    echo "          movie 15 mantle examples"
    echo "----------------------------------------------------------------------------"
    echo "Notes: (1) Run and project per RUN & PROJECT files if unspecified"
    echo "       (2) Variables per VARS file"
    echo "       (3) fps per FPS file if unspecified"
    echo "       (4) Plot  files (*.png) must exist      in $homedir/projects/project/run/plots/"
    echo "       (5) Movie files (*.avi) will be created in $homedir/projects/project/run/movies/"
    echo "       (6) Use the plots alias or run plot_npy.py to create plot files from data files"
    echo "------------------------------------------------------------------------"
    exit
endif

set RUN     = `cat $homedir/RUN`
set PROJECT = `cat $homedir/PROJECT`

if($#argv >= 1) then
    set fps = $1
endif
if($#argv >= 2) then
    set RUN     = $2
endif
if($#argv == 3) then
    set PROJECT = $3
endif
if($#argv == 4) then
    set next_run = $4
else
    set next_run = false
endif

set VARS = `cat $homedir/VARS`

if($next_run == false) then
    echo "============================================================================"
endif
echo "Making movie: run = $RUN, project = "$PROJECT", vars = $VARS"
if($next_run == false) then
    echo "============================================================================"
endif

set nonomatch = 1
set plotdir   = $homedir/projects/$PROJECT/$RUN/plots/
if(!(-e $plotdir)) then
    echo "Plot directory ${plotdir} doesn't exist! Could not make movies: run $RUN, project $PROJECT"
    echo "----------------------------------------------------------------------------"
    exit
endif

set nmovies = 0
foreach what ($VARS)

    set plotfile  = ${plotdir}$what"_0"

    if (-e ${plotfile}000000000*.png) then
        set movfile   = $homedir/projects/$PROJECT/$RUN/movies/$what
        mencoder "mf://${plotfile}*.png" -mf w=800:h=200:fps={$fps}:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o $movfile"_fps="$fps.avi >& /dev/null
        @ nmovies = $nmovies + 1
        if($next_run == false) then
            echo "${movfile}_fps=$fps.avi"
        else
            echo "    ${movfile}_fps=$fps.avi"
        endif
    else
        echo "Plot files do not exist, cannot make movies: var $what, run $RUN, project $PROJECT"
    endif
end
if($next_run == false) then
    echo "============================================================================"
endif
