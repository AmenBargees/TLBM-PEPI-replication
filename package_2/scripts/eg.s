#!/bin/csh -f

if(($#argv == 1 && $1 == "-h") || $#argv > 9) then
    echo "-------------------------------------------------------------------"
    echo "Display using eog up to 9 png file(s)"
    echo "without eog WARNING messages messages"
    echo "-------------------------------------------------------------------"
    echo "Usage   : eg   filename(s)"
    echo "-------------------------------------------------------------------"
    echo "Examples: eg file.png"
    echo "          eg file1.png file2.png"
    echo "-------------------------------------------------------------------"
    echo "Note: eg is an alias to eg.s"
    echo "-------------------------------------------------------------------"
    exit
endif

set i = 1
while($i <= $#argv)
    set file = $argv[$i]
    eog $file >& /dev/null &
    @ i = $i + 1
end
