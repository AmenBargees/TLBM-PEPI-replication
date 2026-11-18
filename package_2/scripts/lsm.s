#!/bin/csh -f

set PROJECT = `cat $homedir/PROJECT`
set RUN     = `cat $homedir/RUN`

if(($#argv == 1 && $1 == "-h") || $#argv > 2) then
    echo "-------------------------------------------------------------------"
    echo 'Do "ls" of the movies directory or a file within it for the run'
    echo "and project specified in the RUN and PROJECT files of your"
    echo "homedir: $homedir"
    echo "-------------------------------------------------------------------"
    echo "Usage:    lsm [option] [filename]"
    echo "          lsm [option] [filename]"
    echo "-------------------------------------------------------------------"
    echo "Examples: lsm                          list movies directory via ls"
    echo '          lsm      u_0000000002_Gy.png list file via ls'
    echo '          lsm -alF T_0000000001_Gy.png list file via ls -alF'
    echo "-------------------------------------------------------------------"
    echo "Note: lsm is an alias to lsm.s"
    echo "-------------------------------------------------------------------"
    exit
endif

set moviedir = $homedir/projects/$PROJECT/$RUN/movies

pushd $moviedir >& /dev/null

if($#argv == 0) then
    ls
else if($#argv == 1) then
    ls $1
else if($#argv == 2) then
    ls $1 $2
endif

popd >& /dev/null
