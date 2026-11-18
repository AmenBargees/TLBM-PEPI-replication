#!/bin/csh -f

if(($#argv == 1 && $1 == -h) || $#argv > 2) then
    echo "-------------------------------------------------------------------"
    echo "more the error file of a run, homedir/projects/project/run.err"
    echo "-------------------------------------------------------------------"
    echo "Usage   :  me [run project]"
    echo "-------------------------------------------------------------------"
    echo "Examples:  me"
    echo "           me archean"
    echo "           me hadean examples"
    echo "-------------------------------------------------------------------"
    echo "Notes: (1) Run & project per RUN and PROJECT files if unspecified"
    echo "       (2) me is an alias to me.s"
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

more $homedir/projects/$PROJECT/$RUN.err
