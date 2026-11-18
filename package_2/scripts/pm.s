#!/bin/csh -f

if(($#argv == 1 && $1 == -h) || $#argv > 2) then
    echo "-------------------------------------------------------------------"
    echo "pushd to movies directory"
    echo "-------------------------------------------------------------------"
    echo "Usage:    pm"
    echo "          pm run"
    echo "          pm run project"
    echo "-------------------------------------------------------------------"
    echo "Examples: pm"
    echo "          pm mantle"
    echo "          pm mantle3D examples"
    echo "-------------------------------------------------------------------"
    echo "Notes: (1) Run & project per RUN and PROJECT files if unspecified"
    echo "       (2) pm is an alias to pm.s"
    echo "-------------------------------------------------------------------"
    exit
endif

set PROJECT = `cat PROJECT`
set RUN     = `cat RUN`

if($#argv > 0) then
    set RUN     = $1
endif
if($#argv > 1) then
    set PROJECT = $2
endif

pushd $homedir/projects/$PROJECT/$RUN/movies/
