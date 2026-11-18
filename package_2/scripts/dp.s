#!/bin/csh -f

set PROJECT = `cat $homedir/PROJECT`
set RUN     = `cat $homedir/RUN`
set plotdir = $homedir/projects/$PROJECT/$RUN/plots

if(($#argv == 1 && $1 == "-h") || $#argv > 9) then
    echo "-------------------------------------------------------------------"
    echo "Display using display up to 9 png file(s) in the plots directory"
    echo "for the run and project specified in the RUN and PROJECT files"
    echo "-------------------------------------------------------------------"
    echo "Usage   : dp   filename(s)"
    echo "-------------------------------------------------------------------"
    echo "Examples: dp T_z_0000000004_Gy.png"
    echo "          dp T_0000000001_Gy.png u_0000000001_Gy.png"
    echo "-------------------------------------------------------------------"
    echo "Note: dp is an alias to dp.s"
    echo "-------------------------------------------------------------------"
    exit
endif

pushd $plotdir >& /dev/null

set i = 1
while($i <= $#argv)
    set file = $argv[$i]
    display $file >& /dev/null &
    @ i = $i + 1
end

popd >& /dev/null
