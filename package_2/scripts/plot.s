#!/bin/csh -f

if(($#argv == 1 && $1 == "-h") || $#argv > 2) then

    echo "-------------------------------------------------------------------"
    echo "Generate plot files from data files for a run of a project"
    echo "-------------------------------------------------------------------"
    echo "Usage:    plot [run project]"
    echo "------------------------------------------------------------------------"
    echo "Examples: plot"
    echo "          plot archean"
    echo "          plot hadean examples"
    echo "-------------------------------------------------------------------"
    echo "Notes: (1) Run and project per RUN & PROJECT files"
    echo "       (2) Data files (*.npy) must exist      in $homedir/projects/project/run/data/"
    echo "       (3) Plot files (*.png) will be created in $homedir/projects/project/run/plots/"
    echo "-------------------------------------------------------------------"
    exit
endif

set PROJECT = `cat $homedir/PROJECT`
set RUN     = `cat $homedir/RUN`

if($#argv == 1) then
    set RUN     = $1
endif
if($#argv == 2) then
    set RUN     = $1
    set PROJECT = $2
endif

pushd $homedir >& /dev/null

echo "--------------------------------------------------------------------"
echo "Making plots from data files of a run of a project, status via any of"
echo "--------------------------------------------------------------------"
echo "tail -20 $homedir/projects/$PROJECT/$RUN/plots/plots.out"
echo "sp"
echo "sp lines $RUN $PROJECT"
echo "--------------------------------------------------------------------"

python3 plot_npy.py $PROJECT $RUN >& ./projects/$PROJECT/$RUN/plots/plots.out &

popd >& /dev/null
