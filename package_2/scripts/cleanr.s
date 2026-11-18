#!/bin/csh -f

if(($#argv == 1 && $1 == "-h") || $#argv > 2) then

    echo "-------------------------------------------------------------------"
    echo "Delete all output files (.out), and data, plot & movie directories"
    echo "for a run of a project. Input files (.in) will not NOT be deleted"
    echo "-------------------------------------------------------------------"
    echo "Usage:    cleanr"
    echo "          cleanr run_name"
    echo "          cleanr run_name project_name"
    echo "-------------------------------------------------------------------"
    echo "Examples: cleanr"
    echo "          cleanr archean"
    echo "          cleanr hadean examples"
    echo "-------------------------------------------------------------------"
    echo "Note: If unspecified, run and project names per RUN & PROJECT files"
    echo "-------------------------------------------------------------------"
    exit
else if($#argv == 0) then
    set RUN     = `cat $homedir/RUN`
    set PROJECT = `cat $homedir/PROJECT`
else if($#argv == 1) then
    set RUN     = $1
    set PROJECT = `cat $homedir/PROJECT`
else if($#argv == 2) then
    set RUN     = $1
    set PROJECT = $2
endif

echo "-------------------------------------------------------------------"
echo "Cleaning $RUN of project $PROJECT"
echo "-------------------------------------------------------------------"

set homedir = $homedir

echo -n "Doing rm $homedir/projects/$PROJECT/$RUN.out ... "
rm $homedir/projects/$PROJECT/$RUN.out        >& /dev/null
echo "Done"

echo -n "Doing rm $homedir/projects/$PROJECT/$RUN.bmk ... "
rm $homedir/projects/$PROJECT/$RUN.bmk        >& /dev/null
echo "Done"

echo -n "Doing rm -r $homedir/projects/$PROJECT/$RUN/plots/ ... "
rm -r $homedir/projects/$PROJECT/$RUN/plots/  >& /dev/null
echo "Done"

echo -n "Doing rm -r $homedir/projects/$PROJECT/$RUN/data/ ... "
rm -r $homedir/projects/$PROJECT/$RUN/data/   >& /dev/null
echo "Done"

echo -n "Doing rm -r $homedir/projects/$PROJECT/$RUN/movies/ ... "
rm -r $homedir/projects/$PROJECT/$RUN/movies/ >& /dev/null
echo "Done"

echo "-------------------------------------------------------------------"
echo "Cleanup complete for $RUN of project $PROJECT"
echo "-------------------------------------------------------------------"
