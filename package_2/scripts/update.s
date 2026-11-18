#!/bin/csh -f

if(($#argv == 1 && $1 == "-h") || $#argv > 1) then
    echo "--------------------------------------------------------------------------"
    echo "Generate GTXLBM_run.py from GTXLBM.py by adding the version number & date,"
    echo "developer & developer, version date, and username so this information will"
    echo "be printed to the output file: homedir/projects/project/run.out"
    echo "--------------------------------------------------------------------------"
    echo "Usage : update"
    echo "--------------------------------------------------------------------------"
    echo "Note  : update is an alias to update.s"
    echo "--------------------------------------------------------------------------"
    exit
endif

set homedir         = $homedir

set release_author    = `cat $homedir/.author`
set release_version   = `cat $homedir/release_version`
set release_date      = `cat $homedir/release_date`

set version_developer = `cat $homedir/version_developer`
set version_date      = `cat $homedir/version_date`
set version_user      = `whoami`

set todays_date     = `$homedir/scripts/get_date.s`

echo "=========================================================================="
if(-e $homedir/.author && $version_user == peter) then
    set version_developer = $release_author
    set update    = true
    if($#argv == 1) then

        set  release_version    = $1
        set  release_date       = $todays_date

        set  version_developer  = $release_author
        set  version_date       = $todays_date

        echo "New release of GTXLBM created with the following information"
        echo "=========================================================================="
        echo "release_author      : $release_author"
        echo "release_version     : $release_version"
        echo "release_date        : $release_date"
        echo "--------------------------------------------------------------------------"
        echo "version_developer   : $version_developer"
        echo "version_date        : $version_date"
    endif
else
    set version_developer = $version_user
    set version_date      = $todays_date

    echo "A new version of GTXLBM will be created with the following information"
    echo "=========================================================================="
    echo "release_author    : $release_author"
    echo "release_version   : $release_version"
    echo "release_date      : $release_date"
    echo "--------------------------------------------------------------------------"
    echo "version_developer : $version_developer"
    echo "version_date      : $version_date"
    echo "--------------------------------------------------------------------------"
    echo -n "Do you want to contnue (y or n)? "
    set response = $<
    if($response == y) then
        set update = true
    else
        set update = false
    endif

    echo "--------------------------------------------------------------------------"
    echo "New version of GTXLBM created with the following information"
    echo "=========================================================================="
    echo "release_author    : $release_author"
    echo "release_version   : $release_version"
    echo "release_date      : $release_date"
    echo "--------------------------------------------------------------------------"
    echo "version_developer : $version_developer"
    echo "version_date      : $version_date"
endif

pushd $homedir >& /dev/null

echo $release_author    > $homedir/release_author
echo $release_version   > $homedir/release_version
echo $release_date      > $homedir/release_date

echo $version_developer > $homedir/version_developer
echo $version_date      > $homedir/version_developer

if($update == true) then

    sed s/RELEASE_VERSION/$release_version/ $homedir/GTXLBM.py | sed s/RELEASE_DATE/$release_date/ | sed s/RELEASE_AUTHOR/$release_author/ | sed s/VERSION_DEVELOPER/$version_developer/ | sed s/VERSION_USER/$version_user/ | sed s/VERSION_DATE/$version_date/ > $homedir/GTXLBM_run.py
    echo "GTXLBM_run.py was created with the following information"
    echo "========================================================================"
    echo "release_author    : $release_author"
    echo "release_version   : $release_version"
    echo "release_date      : $release_date"
    echo "--------------------------------------------------------------------------"
    echo "version_developer : $version_developer"
    echo "version_date      : $version_date"
    echo "========================================================================"
    echo "GTXLBM_run.py is ready to run"
else
    echo "GTXLBM_run.py has not been updated"
endif
echo "========================================================================"

popd >& /dev/null
