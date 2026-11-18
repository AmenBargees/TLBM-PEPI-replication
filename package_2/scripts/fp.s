#!/bin/csh -f

set PROJECT = `cat $homedir/PROJECT`
set RUN     = `cat $homedir/RUN`
set plotdir = $homedir/projects/$PROJECT/$RUN/plots

if(($#argv == 1 && $1 == "-h") || $#argv > 9) then
    echo "-------------------------------------------------------------------"
    echo "Display using fbi up to 9 png file(s) in the plots directory"
    echo "for the run and project specified in the RUN and PROJECT"
    echo "files of your homedir: $homedir"
    echo "-------------------------------------------------------------------"
    echo "Usage   : fp   filename(s)"
    echo "-------------------------------------------------------------------"
    echo "Examples: fp T_z_0000000004_Gy.png"
    echo "          fp T_0000000001_Gy.png u_0000000001_Gy.png"
    echo "-------------------------------------------------------------------"
    echo "Note: fp is an alias to fp.s"
    echo "-------------------------------------------------------------------"
    exit
endif

pushd $plotdir >& /dev/null

set i = 1
while($i <= $#argv)
    set file = $argv[$i]
    fbi $file &

    @ i = $i + 1
end

popd >& /dev/null
