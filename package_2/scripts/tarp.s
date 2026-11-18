#!/bin/csh -f

set homedir         = $homedir
set user            = `whoami`
set release_version = `cat $homedir/release_version`
set release_date    = `cat $homedir/release_date`
set release_author  = `cat $homedir/.author`
set release         = "GTXLBM_V${release_version}_${release_date}"

if(($#argv == 1 && $1 == "-h") || ($#argv == 2 && $1 != "time") || $#argv > 2) then
    echo "------------------------------------------------------------------------"
    if(($#argv == 2 && $1 != "time") || $#argv > 2) then
        echo "------------------------------------------------------------------------"
        echo "ERROR: Incorrect usage of tarp, see tarp documentation below"
    endif
    echo "------------------------------------------------------------------------"
    echo "Create tarfile named as below of GTXLBM source code, scripts, and a project"
    echo "-------------------------------------------------------------------"
    echo "homedir/tarfiles/GTXLBM_V{release_version}_{release_date}_{user}_{date}.tar"
    echo "-------------------------------------------------------------------"
    echo "Usage:    tarp [arg1,[arg2]]"
    echo ""
    echo "If only one argument  is  specified, arg1 must be time          or      a project name"
    echo "If two      arguments are specified, arg1 must be time and arg2 must be a project name"
    echo "If the project name is unspecified, it will be obtained from the PROJECT file"
    echo "-------------------------------------------------------------------"
    echo "Examples:"
    echo "    tarp               (tarfile of a project with today's date        where project is per PROJECT file)"
    echo "    tarp project       (tarfile of a project with today's date)"
    echo "    tarp time          (tarfile of a project with today's date & time where project is per PROJECT file)"
    echo "    tarp time project  (tarfile of a project with today's date & time)"
    echo "------------------------------------------------------------------------"
    echo "release_version & release_date are the GTXLBM version & its release date"
    echo "------------------------------------------------------------------------"
    echo "homedir = $homedir is where you installed GTXLBM"
    echo "user    = $user is your user name per whoami"
    echo "date    = YYYY-MM-DD or YYYY-MM-DD-HH:MM if the time option is specified"
    echo "project = as specified, or if unspecified, via the PROJECT file"
    echo "------------------------------------------------------------------------"
    echo "Note: tarp is an alias to tarp.s"
    echo "------------------------------------------------------------------------"
    exit
endif

set project = `cat ${homedir}/PROJECT`

if($#argv == 1 && $1 == "time") then
    set date = `$homedir/scripts/get_date.s time`
else
    set date = `$homedir/scripts/get_date.s`
    echo "nargv = $#argv"
    if($#argv == 1) then
        set project = $1
        echo "Set project to arg1 = $1"
    else if($#argv == 2) then
        set project = $2
        echo "Set project to arg2 = $2"
    endif
endif

set n_saved_projects = 0

if(-e $homedir/.author && $user == peter) then
    set version = "GTXLBM_V${release_version}_${release_author}_${date}"
else
    set version = "GTXLBM_V${release_version}_${user}_${date}"
endif

set tarfile = $homedir/tarfiles/$version.tar

pushd $homedir >& /dev/null

mkdir -p tarfiles

# Initialize the number of saved and missing files and/or directories

set files = `cat .files`

echo "========================================================================"
echo "Creating ${tarfile}"
echo "========================================================================"
echo "Source code    : *.py, *.c"
echo "------------------------------------------------------------------------"
echo "Files:         : $files"
echo "------------------------------------------------------------------------"
echo "papers         : `ls papers`"
echo "------------------------------------------------------------------------"
echo "C-shell scripts: ./scripts/*.s"
echo "------------------------------------------------------------------------"
echo "project        : $project"
echo "========================================================================"

echo -n "Adding source code and shared object file ... "
tar cvf $tarfile *.py *.c *.so .author --warning=no-file-changed >& /dev/null
echo "Done!"

echo -n "Adding files ... "
tar rf $tarfile $files --warning=no-file-changed >& /dev/null
echo "Done!"

echo -n "Adding papers ... "
tar rf $tarfile papers --warning=no-file-changed
echo "Done!"

echo -n "Adding scripts ... "
tar rf $tarfile scripts --warning=no-file-changed
echo "Done!"

set project_dir = ./projects/$project/
echo -n "Adding project $project ... "
if (!(-e $project_dir)) then
    echo "Failed! Project directory does not exist."
else
    @ n_saved_projects = $n_saved_projects + 1
    tar rf $tarfile  $project_dir --warning=no-file-changed
    echo "Done!"
endif

if($n_saved_projects == 0) then
    echo "========================================================================"
    echo "ERROR: The project directory not exist so it was not saved!"
    echo "========================================================================"
    echo "Missing project directory: $project_dir"
    echo "Check project name was specified correctly and its directory exists."
else
    echo "========================================================================"
    echo "Tarfile created: "`du -s -h ${tarfile}`
endif
echo "========================================================================"

popd >& /dev/null
