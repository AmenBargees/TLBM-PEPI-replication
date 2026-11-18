#!/bin/csh -f

if(($#argv == 1 && $1 == "-h") || ($#argv == 1 && $1 != "time") || $#argv > 1) then
    echo "-------------------------------------------------------------------"
    echo "Get today's date or today's date and time"
    echo "-------------------------------------------------------------------"
    echo "Usage:    get_date [time]"
    echo "-------------------------------------------------------------------"
    echo "Examples: get_date       Get today's date          in format YYYY-MM-DD"
    echo "          get_date time  Get today's date and time in format YYYY-MM-DD-HH:MM"
    echo "-------------------------------------------------------------------"
    echo "Note: get_date is an alias to get_date.s"
    echo "-------------------------------------------------------------------"
    exit
endif

if($1 == "time") then
    set date = `date +"%F-%H:%M"`
else
    set date = `date -Idate`
endif

echo $date

