#!/bin/csh -f

if(($#argv == 1 && $1 == -h) || $#argv > 2) then
    echo "-------------------------------------------------------------------"
    echo "more the output file of a plot generation run:"
    echo "    homedir/projects/project/run/plots/plots.out"
    echo "-------------------------------------------------------------------"
    echo "Usage   :  mp [run project]"
    echo "-------------------------------------------------------------------"
    echo "Examples:  mp"
    echo "           mp archean"
    echo "           mp mantle examples"
    echo "-------------------------------------------------------------------"
    echo "Notes: (1) Run & project per RUN and PROJECT files if unspecified"
    echo "       (2) mp is an alias to mp.s" 
    echo "-------------------------------------------------------------------"
    exit
endif

set RUN     = `cat $homedir/RUN`
set PROJECT = `cat $homedir/PROJECT`

if($#argv > 0) then
    set RUN     = $1
endif

if($#argv > 1) then
    set PROJECT = $2
endif

more $homedir/projects/$PROJECT/$RUN/plots/plots.out
