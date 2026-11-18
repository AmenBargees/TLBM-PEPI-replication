#!/bin/csh -f

if(($#argv == 1 && $1 == -h) || $#argv > 2) then
    echo "-------------------------------------------------------------------"
    echo "more the output file of a run, homedir/projects/project/run.out"
    echo "-------------------------------------------------------------------"
    echo "Usage   :  mo [run project]"
    echo "-------------------------------------------------------------------"
    echo "Examples:  mo"
    echo "           mo archean"
    echo "           mo hadean examples"
    echo "-------------------------------------------------------------------"
    echo "Notes: (1) Run & project per RUN and PROJECT files if unspecified"
    echo "       (2) mo is an alias to mo.s"
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

more $homedir/projects/$PROJECT/$RUN.out
