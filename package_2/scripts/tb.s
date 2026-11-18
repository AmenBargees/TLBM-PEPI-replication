#!/bin/csh -f

set lines = 200

if(($#argv == 1 && $1 == "-h") || $#argv > 3) then
    echo "-------------------------------------------------------------------"
    echo "Tail the run.out file of a GTXLBM run"
    echo "    homedir/projects/project/run.out"
    echo "-------------------------------------------------------------------"
    echo "Usage:    tb [lines run project]"
    echo "--------------------------------------------------------------------------"
    echo "Examples: tb"
    echo "          tb 100"
    echo "          tb 50 archean"
    echo "          tb 40 hadean examples"
    echo "          tb 80 mantle examples"
    echo "--------------------------------------------------------------------------"
    echo "Notes: (1) Run & project per RUN and PROJECT files if unspecified"
    echo "       (2) Lines defaults tb $lines if unspecified"
    echo "       (3) tb is an alias tb tb.s"
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

tail -$lines $homedir/projects/$PROJECT/$RUN.bmk
