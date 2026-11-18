#!/bin/csh -f

set PROJECT = `cat $homedir/PROJECT`
set RUN     = `cat $homedir/RUN`

if(($#argv == 1 && $1 == "-h") || $#argv > 2) then
    echo "-------------------------------------------------------------------"
    echo 'Do "ls" of the plots directory or a file within it for the run'
    echo "and project specified in the RUN and PROJECT files of your"
    echo "homedir: $homedir"
    echo "-------------------------------------------------------------------"
    echo "Usage:    lsp   [option] [filename]"
    echo "          lsp   name \*"
    echo "          lsp   \* \name"
    echo "-------------------------------------------------------------------"
    echo "Examples: lsp                          list plots directory via ls"
    echo '          lsp      u_0000000100_Gy.png list file'
    echo '          lsp -alF T_0000000100_Gy.png list file via ls -alF'
    echo '          lsp T_z \*                   list files matching T_z*'
    echo '          lsp \* 0_Gy.png              list files matching *0_Gy.png'
    echo '          lsp \* 00_My.png             list files matching *00_My.png'
    echo "-------------------------------------------------------------------"
    echo "Note: lsp is an alias to lsp.s"
    echo "-------------------------------------------------------------------"
    exit
endif

set plotdir = $homedir/projects/$PROJECT/$RUN/plots

pushd $plotdir >& /dev/null
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
