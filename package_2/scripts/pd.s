#!/bin/csh -f

if(($#argv == 1 && $1 == -h) || $#argv > 2) then
    echo "-------------------------------------------------------------------"
    echo "pushd to data directory"
    echo "-------------------------------------------------------------------"
    echo "Usage:    pd"
    echo "          pd run"
    echo "          pd run project"
    echo "-------------------------------------------------------------------"
    echo "Examples: pd"
    echo "          pd mantle"
    echo "          pd mantle3D examples"
    echo "-------------------------------------------------------------------"
    echo "Notes: (1) Run & project per RUN and PROJECT files if unspecified"
    echo "       (2) pd is an alias to pd.s"
    echo "-------------------------------------------------------------------"
    exit
endif

set PROJECT = `cat $homedir/PROJECT`
set RUN     = `cat $homedir/RUN`

if($#argv > 0) then
    set RUN     = $1
endif
if($#argv > 1) then
    set PROJECT = $2
endif

pushd $homedir/projects/$PROJECT/$RUN/data/
