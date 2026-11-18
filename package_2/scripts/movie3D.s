#!/bin/csh -f

set fps = `cat $homedir/FPS`

if($#argv < 2 || $#argv > 5) then
    echo "----------------------------------------------------------------------------"
    echo "Make movies from plot files of a run of a project from a 3D GTXLBM run"
    echo "----------------------------------------------------------------------------"
    echo "Usage:    movie3D.s plane coordinates [fps run project]"
    echo "----------------------------------------------------------------------------"
    echo "Examples: movie3D.s xy 10,50"
    echo "          movie3D.s xz 20,50 20"
    echo "          movie3D.s xy 6     10 mantle3D"
    echo "          movie3D.s xy 6,54  15 mantle3D examples"
    echo "----------------------------------------------------------------------------"
    echo "Notes: (1) Run and project per RUN & PROJECT files if unspecified"
    echo "       (2) Variables per VARS file"
    echo "       (3) fps per FPS file if unspecified"
    echo "       (4) Plot  files (*.png) must exist      in $homedir/projects/project/run/plots/"
    echo "       (5) Movie files (*.avi) will be created in $homedir/projects/project/run/movies/"
    echo "       (6) Use the plots alias or run plot_npy.py to create plot files from data files"
    echo "----------------------------------------------------------------------------"
    exit
endif

set RUN         = `cat $homedir/RUN`
set PROJECT     = `cat $homedir/PROJECT`
set PLANE       = $1
set COORDINATES = ($2:as/,/ /)

if($#argv > 2) then
    set fps     = $3
endif
if($#argv > 3) then
    set RUN     = $4
endif
if($#argv > 4) then
    set PROJECT = $5
endif

if($PLANE == xy) then
    set COORD = z 
else if ($PLANE == xz) then
    set COORD = y 
else if ($PLANE == xy) then
    set COORD = z
endif

set VARS = `cat $homedir/VARS`

echo "============================================================================"
echo "Making movies: PROJECT = "$PROJECT", RUN = "$RUN" at "$fps" frames per second"
echo "plane = ${PLANE}, coordinates = $2, variables per VARS file = "$VARS
echo "----------------------------------------------------------------------------"

set nonomatch = 1
set plotdir   = $homedir/projects/$PROJECT/$RUN/plots/

if(-e $plotdir) then
    echo "plotdir = ${plotdir} exists, making movies ..."
    echo "----------------------------------------------------------------------------"
else
    echo "plotdir = ${plotdir} doesn't exist! No movies can be made."
    echo "----------------------------------------------------------------------------"
    exit
endif

set VARS = `cat $homedir/VARS`

set nmovies = 0
foreach what ($VARS)

    foreach coordinate ($COORDINATES)
        set SPECS       = "_plane=${PLANE}_${COORD}=${coordinate}"

        set plotfile  = ${plotdir}$what$SPECS"_"

        echo -n "Checking plot files exist for "$what" ... "
        set plotfile0 = ${plotfile}0000000000_*.png || true

        if(-e ${plotfile0}) then
            set movfile   = $homedir/projects/$PROJECT/$RUN/movies/$what$SPECS
            echo -n "making movie file ${movfile}_fps=$fps.avi ... "
            mencoder "mf://${plotfile}*.png" -mf w=800:h=200:fps={$fps}:type=png -ovc lavc -lavcopts vcodec=mpeg4:mbd=2:trell -oac copy -o $movfile"_fps="$fps.avi >& /dev/null
            echo "Done"
            @ nmovies = $nmovies + 1
        else
            echo "plot files do not exist, skipping ${what}"
        endif
    end
end
echo "----------------------------------------------------------------------------"
if($nmovies == 0) then
    echo "There are no plot files from which to make movies for this project and run"
else
    echo "Finished making all "$nmovies" movies"
endif
echo "----------------------------------------------------------------------------"
