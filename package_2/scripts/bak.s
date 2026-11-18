#!/bin/csh -f

set homedir         = $homedir
set user            = `whoami`
set release_version = `cat $homedir/release_version`
set release_date    = `cat $homedir/release_date`
set release_author  = `cat $homedir/.author`
set release         = "GTXLBM_V${release_version}_${release_date}"

if(($#argv == 1 && $1 == "-h") || ($#argv == 1 && $1 != "time") || ($#argv > 1)) then
    echo "------------------------------------------------------------------------"
    echo "Create backup directory of GTXLBM source code, scripts, and the projects"
    echo "and runs specified in the PROJECTS, RUNS & DIRS files per formats below"
    echo "-------------------------------------------------------------------"
    echo "Usage:    bak.s [time]"
    echo "-------------------------------------------------------------------"
    echo "Examples: bak.s       (Create backup directory with today's date)"
    echo "          bak.s time  (Create backup directory with today's date & time)"
    echo "-------------------------------------------------------------------"
    echo "Backup directory name:"
    echo "${homedir}/backups/GTXLBM_V{release_version}_{release_date}_{user}_{date}"
    echo "------------------------------------------------------------------------"
    echo "Date format: YYYY-MM-DD or YYYY-MM-DD-HH:MM if the time option is specified"
    echo "------------------------------------------------------------------------"
    echo "User per whoami, release_version & release_date are GTXLBM version & its release date"
    echo "------------------------------------------------------------------------"
    echo "Filenames & formats specifying backup directory contents are as below"
    echo "------------------------------------------------------------------------"
    echo "PROJECTS file to specify project(s), one line as below"
    echo ""
    echo "project1 project2 ..."
    echo ""
    echo "RUNS file to specify run(s) from project(s), one line per project as below"
    echo ""
    echo "project1 run1_of_project1 run2_of_project2 ..."
    echo "project2 run1_of_project2 run2_of_project2 ..."
    echo ""
    echo "DIRS file to specify plots and/or data and/or bmk directories from the"
    echo "specified project(s) and run(s), one line as below"
    echo ""
    echo "plots movies data bmk"
    echo "------------------------------------------------------------------------"
    exit
endif

if($#argv == 1 && $1 == "time") then
    set date = `$homedir/scripts/get_date.s time`
else
    set date = `$homedir/scripts/get_date.s`
endif

if(-e $homedir/.author && $user == peter) then
    set version = "GTXLBM_V${release_version}_${release_author}_${date}"
else
    set version = "GTXLBM_V${release_version}_${user}_${date}"
endif

set bakdir  = $homedir/backups/$version
set tarfile = $homedir/backups/$version.tar

if($#argv == 1 && $1 == "time") then
    tar.s bak time
else
    tar.s bak
endif

echo "Creating backup directory: $bakdir"

mkdir -p $bakdir

if(-e $bakdir && !(-e $homedir/.author && $user == peter)) then
    echo "========================================================================"
    echo "WARNING: $bakdir"
    echo ""
    echo "         already exists and will be overwritten if you continue"
    echo "------------------------------------------------------------------------"
    echo -n "Do you want to continue (y or n)? "
    set response = $<
    if($response == y) then
        echo "Removing old backup dir  : $bakdir"
        rm -r $bakdir
    else
        echo ""
        echo "Exiting!"
        echo ""
        echo "Backup dir unchanged     : $bakdir"
        echo "========================================================================"
        exit
    endif
endif

pushd $bakdir >& /dev/null

echo "ls -alF tarfile: "`ls -alF $tarfile`

echo "------------------------------------------------------------------------"
echo "Extracting tarfile       : $tarfile"
echo "------------------------------------------------------------------------"

tar xvf $tarfile >& /dev/null

echo "Removing tarfile: $tarfile"

rm $tarfile

echo "========================================================================"
echo "Backup directory created : `du -s -h $bakdir`"
echo "========================================================================"
popd >& /dev/null
