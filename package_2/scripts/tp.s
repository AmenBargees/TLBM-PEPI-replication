#!/bin/csh -f

set lines = 15

if(($#argv == 1 && $1 == "-h") || $#argv > 3) then
    echo "-------------------------------------------------------------------"
    echo "Tail the plots.out file of a plot_npy.py run"
    echo "    homedir/projects/project/run/plots/plots.out"
    echo "-------------------------------------------------------------------"
    echo "Usage:    tp [lines run project]"
    echo "--------------------------------------------------------------------------"
    echo "Examples: tp"
    echo "          tp 30"
    echo "          tp 100 archean"
    echo "          tp 20 hadean examples"
    echo "--------------------------------------------------------------------------"
    echo "Notes: (1) Run & project per RUN and PROJECT files if unspecified"
    echo "       (2) Lines defaults to $lines if unspecified"
    echo "       (3) tp is an alias to tp.s"
    echo "-------------------------------------------------------------------"
    exit
endif

set homedir = $homedir
set RUN     = `cat $homedir/RUN`
set PROJECT = `cat $homedir/PROJECT`

if($#argv >= 1) then
    set lines = $1
endif
if($#argv >= 2) then
    set RUN     = $2
endif
if($#argv == 3) then
    set PROJECT = $3
endif

tail -$lines $homedir/projects/$PROJECT/$RUN/plots/plots.out
