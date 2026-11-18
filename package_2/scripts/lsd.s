#!/bin/csh -f

set PROJECT = `cat $homedir/PROJECT`
set RUN     = `cat $homedir/RUN`

if(($#argv == 1 && $1 == "-h") || $#argv > 2) then
    echo "-------------------------------------------------------------------"
    echo 'Do "ls" of the data directory or a file within it for the run'
    echo "and project specified in the RUN and PROJECT files of your"
    echo "your homedir: $homedir"
    echo "-------------------------------------------------------------------"
    echo "Usage   : lsd   [ls option] [filename]"
    echo "-------------------------------------------------------------------"
    echo "Examples: lsd                          list data directory via ls"
    echo '          lsd      u_0000000002_Gy.png list file via ls'
    echo '          lsd -alF T_0000000001_Gy.png list file via ls -alF'
    echo "-------------------------------------------------------------------"
    echo "Note: lsd is an alias to lsd.s"
    echo 'variable set by "source homedir/cshrc_additions"'
    echo "-------------------------------------------------------------------"
    exit
endif

set datadir = $homedir/projects/$PROJECT/$RUN/data

pushd $datadir >& /dev/null

if($#argv == 0) then
    ls
else if($#argv == 1) then
    ls $1
else if($#argv == 2) then
    ls $1 $2
endif

popd >& /dev/null
