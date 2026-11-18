#!/bin/csh -f

if(($#argv == 1 && $1 == -h) || $#argv > 1) then
    echo "-------------------------------------------------------------------"
    echo "pushd to bmk directory of a project"
    echo "-------------------------------------------------------------------"
    echo "Usage:    pb"
    echo "          pb project"
    echo "-------------------------------------------------------------------"
    echo "Examples: pb"
    echo "          pb bencharks_annulus"
    echo "-------------------------------------------------------------------"
    echo "Notes: (1) project per PROJECT file if unspecified"
    echo "       (2) pb is an alias to pb.s"
    echo "-------------------------------------------------------------------"
    exit
endif

set PROJECT = `cat $homedir/PROJECT`

if($#argv > 0) then
    set PROJECT = $1
endif

pushd $homedir/projects/$PROJECT/bmk/
