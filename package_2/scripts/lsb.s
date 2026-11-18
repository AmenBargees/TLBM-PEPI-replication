#!/bin/csh -f

set PROJECT = `cat $homedir/PROJECT`

if(($#argv == 1 && $1 == "-h") || $#argv > 2) then
    echo "-------------------------------------------------------------------"
    echo 'Do "ls" of the bmk directory or a file within it for the project'
    echo "specified in the PROJECT file of your homedir: $homedir"
    echo "-------------------------------------------------------------------"
    echo "Usage:    lsb   [option] [filename]"
    echo "          lsb   name \*"
    echo "          lsb   \* \name"
    echo "-------------------------------------------------------------------"
    echo "Examples: lsb                          list bmk directory via ls"
    echo '          lsb filename.png             list file named filename.png'
    echo '          lsb -alF filename.png        list file via ls -alF'
    echo '          lsb Pr=1000 \*               list files matching Pr=1000*'
    echo '          lsb \* Nu.png                list files matching *Nu.png'
    echo "-------------------------------------------------------------------"
    echo "Note: lsb is an alias to lsb.s"
    echo "-------------------------------------------------------------------"
    exit
endif

set bmkdir = $homedir/projects/$PROJECT/bmk

pushd $bmkdir >& /dev/null
set noglob

if($#argv == 0) then
    unset noglob
    ls
else if($#argv == 1) then
    unset noglob
    ls $1
else if($#argv == 2 && "$1" == "*") then
    unset noglob
    ls *"$2"
else if($#argv == 2 && "$2" == "*") then
    unset noglob
    ls "$1"*
else if($#argv == 2) then
    unset noglob
    ls $1 $2
endif

popd >& /dev/null
