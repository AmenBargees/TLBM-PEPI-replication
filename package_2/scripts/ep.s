#!/bin/csh -f

set PROJECT = `cat $homedir/PROJECT`
set RUN     = `cat $homedir/RUN`
set plotdir = $homedir/projects/$PROJECT/$RUN/plots

if(($#argv == 1 && $1 == "-h") || $#argv > 9) then
    echo "-------------------------------------------------------------------"
    echo "Display using eog up to 9 png file(s) in the plots directory"
    echo "for the run and project specified in the RUN and PROJECT"
    echo "files of your homedir: $homedir"
    echo "-------------------------------------------------------------------"
    echo "Usage   : ep   filename(s)"
    echo "-------------------------------------------------------------------"
    echo "Examples: ep T_z_0000000004_Gy.png"
    echo "          ep T_0000000001_Gy.png u_0000000001_Gy.png"
    echo "-------------------------------------------------------------------"
    echo "Note: ep is an alias to ep.s"
    echo "-------------------------------------------------------------------"
    exit
endif

pushd $plotdir >& /dev/null

set i = 1
while($i <= $#argv)
    set file = $argv[$i]
    eog $file >& /dev/null &
    @ i = $i + 1
end

popd >& /dev/null
