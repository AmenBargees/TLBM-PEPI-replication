#!/bin/csh -f

set lines = 20

if(($#argv == 1 && $1 == "-h") || $#argv > 3) then
    echo "-------------------------------------------------------------------"
    echo "Tail the standard error file of GTXLBM run done via run or runs"
    echo "    home/projects/project/run.err"
    echo "-------------------------------------------------------------------"
    echo "Usage:    te [lines run project]"
    echo "--------------------------------------------------------------------------"
    echo "Examples: te"
    echo "          te 100"
    echo "          te 50 archean"
    echo "          te 40 hadean examples"
    echo "          te 80 mantle examples"
    echo "--------------------------------------------------------------------------"
    echo "Notes: (1) Run & project per RUN and PROJECT files if unspecified"
    echo "       (2) Lines defaults to $lines if unspecified"
    echo "       (3) te is an alias to te.s"
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

tail -$lines $homedir/projects/$PROJECT/$RUN.err
