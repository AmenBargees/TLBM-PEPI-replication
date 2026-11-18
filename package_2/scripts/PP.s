#!/bin/csh -f

if(($#argv == 1 && $1 == -h) || $#argv > 1) then
    echo "-------------------------------------------------------------------"
    echo "pushd to project directory"
    echo "-------------------------------------------------------------------"
    echo "Usage: pproject  [project]"
    echo "-------------------------------------------------------------------"
    echo "Examples: pproj"
    echo "          pproj examples"
    echo "-------------------------------------------------------------------"
    echo "Notes: (1) Project per PROJECT file if unspecified"
    echo "       (2) PP is an alias to PP.s"
    echo "-------------------------------------------------------------------"
    exit
endif

if($#argv == 0) then
    set PROJECT = `cat $homedir/PROJECT`
else if($#argv == 1) then
    set PROJECT = $1
endif

pushd $homedir/projects/$PROJECT/
