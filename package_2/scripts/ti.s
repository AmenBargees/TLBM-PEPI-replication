#!/bin/csh -f

set lines = 200

if(($#argv == 1 && $1 == "-h") || $#argv > 3) then
    echo "-------------------------------------------------------------------"
    echo "Tail the run.in file of a GTXLBM run"
    echo "    homedir/projects/project/run.in"
    echo "-------------------------------------------------------------------"
    echo "Usage:    ti [lines run project]"
    echo "--------------------------------------------------------------------------"
    echo "Examples: ti"
    echo "          ti 100"
    echo "          ti 50 archean"
    echo "          ti 40 hadean examples"
    echo "          ti 80 mantle examples"
    echo "--------------------------------------------------------------------------"
    echo "Notes: (1) Run & project per RUN and PROJECT files if unspecified"
    echo "       (2) Lines defaults to $lines if unspecified"
    echo "       (3) ti is an alias to ti.s"
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

tail -$lines $homedir/projects/$PROJECT/$RUN.in
