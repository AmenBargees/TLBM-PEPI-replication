#!/bin/csh -f

set lines = 200

if(($#argv == 1 && $1 == "-h") || $#argv > 3) then
    echo "-------------------------------------------------------------------"
    echo "Tail the run.out file of a GTXLBM run"
    echo "    homedir/projects/project/run.out"
    echo "-------------------------------------------------------------------"
    echo "Usage:    to [lines run project]"
    echo "--------------------------------------------------------------------------"
    echo "Examples: to"
    echo "          to 100"
    echo "          to 50 archean"
    echo "          to 40 hadean examples"
    echo "          to 80 mantle examples"
    echo "--------------------------------------------------------------------------"
    echo "Notes: (1) Run & project per RUN and PROJECT files if unspecified"
    echo "       (2) Lines defaults to $lines if unspecified"
    echo "       (3) to is an alias to to.s"
    echo "-------------------------------------------------------------------"
    exit
endif

set RUN     = `cat $homedir/RUN`
set PROJECT = `cat $homedir/PROJECT`

if($#argv >= 1) then
    set lines   = $1
endif
if($#argv >= 2) then
    set RUN     = $2
endif
if($#argv == 3) then
    set PROJECT = $3
endif

tail -$lines $homedir/projects/$PROJECT/$RUN.out
