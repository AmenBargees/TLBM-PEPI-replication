#!/bin/csh -f

set PROJECT = `cat $homedir/PROJECT`
set RUN     = `cat $homedir/RUN`
set plotdir = $homedir/projects/$PROJECT/$RUN/plots

if(($#argv == 1 && $1 == "-h") || $#argv > 9) then
    echo "-------------------------------------------------------------------"
    echo "Display using catimg up to 9 png file(s) in the plots directory"
    echo "for the run and project specified in the RUN and PROJECT files"
    echo "of your homedir: $homedir"
    echo "-------------------------------------------------------------------"
    echo "Usage   : Cp filename(s)"
    echo "-------------------------------------------------------------------"
    echo "Examples: Cp T_z_0000000004_Gy.png"
    echo "          Cp T_0000000001_Gy.png u_0000000001_Gy.png"
    echo "-------------------------------------------------------------------"
    echo "Note: Cp is an alias to Cp.s"
    echo "and it assumes homedir/scripts was added to your PATH envoronmental"
    echo 'variable set by "source homedir/cshrc_additions"'
    echo "-------------------------------------------------------------------"
    echo "WARNING: catimg cat's a low resolution image to your terminal window"
    echo "so you will lose this window. If this happens, open another terminal"
    echo "-------------------------------------------------------------------"
    echo "NOT RECOMMENDED as you will lose your terminal window"
    echo "-------------------------------------------------------------------"
    exit
endif

pushd $plotdir >& /dev/null

set i = 1
while($i <= $#argv)
    set file = $argv[$i]
    catimg $file >& /dev/null &
    @ i = $i + 1
end

popd >& /dev/null
