#!/bin/csh -f

set homedir         = $homedir
set user            = `whoami`
set release_version = `cat $homedir/release_version`
set release_date    = `cat $homedir/release_date`
set release_author  = `cat $homedir/.author`
set release         = "GTXLBM_V${release_version}_${release_date}"

if(($#argv == 1 && $1 == "-h") || ($#argv == 1 && ($1 != "bak" && $1 != "time")) || ($#argv == 2 && $2 != "time") || ($#argv > 4)) then
    echo "------------------------------------------------------------------------"
    echo "Create tarfile named as below of GTXLBM source code, scripts, and: (1) all"
    echo "runs/projects, (2) a specified run of a project, or a specified project"
    echo "-------------------------------------------------------------------"
    echo "homedir/tarfiles/GTXLBM_V{release_version}_{release_date}_{user}_{date}.tar"
    echo "-------------------------------------------------------------------"
    echo "Usage:    tar.s [arg1,[arg2]]"
    echo "-------------------------------------------------------------------"
    echo "Examples:"
    echo "    tar.s all       (tarfile of all runs/projects with today's date)"
    echo "    tar.s all time  (tarfile of all runs/projects with today's date & time )"
    echo "    tar.s run       (tarfile of run  per RUN  & PROJECT  files with today's date)"
    echo "    tar.s run time  (tarfile of run  per RUN  & PROJECT  files with today's date & time )"
    echo "    tar.s runs      (tarfile of runs per RUNS & PROJECTS files with today's date)"
    echo "    tar.s runs time (tarfile of runs per RUNS & PROJECTS files with today's date & time )"
    echo "------------------------------------------------------------------------"
    echo "Where homedir = $homedir"
    echo "Date format   = YYYY-MM-DD or YYYY-MM-DD-HH:MM if the time option is specified"
    echo "------------------------------------------------------------------------"
    echo "User per whoami, release_version & release_date are the GTXLBM version & its release date"
    echo "------------------------------------------------------------------------"
    echo "Filenames/formats specifying tarfile contents if run or runs is specified:"
    echo "------------------------------------------------------------------------"
    echo "PROJECT(S) file to specify project(s), one line as below"
    echo ""
    echo "project1 project2 ..."
    echo ""
    echo "RUN(S) file to specify run(s) from project(s), one line per project as below"
    echo ""
    echo "project1 run1_of_project1 run2_of_project2 ..."
    echo "project2 run1_of_project2 run2_of_project2 ..."
    echo ""
    echo "DIRS file to specify directories (ie. plots and/or data and/or bmk directories)"
    echo "from the specified project(s) and run(s), one line such as as below"
    echo ""
    echo "plots movies data bmk"
    echo "------------------------------------------------------------------------"
    exit
endif

if($1 == "bak") then
    set bak  = true
    if($#argv == 2 && $2 == "time") then
        set date = `$homedir/scripts/get_date.s time`
    else
        set date = `$homedir/scripts/get_date.s`
    endif
else
    set bak  = false
    if($#argv == 1 && $1 == "time") then
        set date = `$homedir/scripts/get_date.s time`
    else
        set date = `$homedir/scripts/get_date.s`
    endif
endif

if(-e $homedir/.author && $user == peter) then
    set version = "GTXLBM_V${release_version}_${release_author}_${date}"
else
    set version = "GTXLBM_V${release_version}_${user}_${date}"
endif

pushd $homedir >& /dev/null

if($bak == true) then
    mkdir -p backups
    set tarfile = $homedir/backups/${version}.tar
else
    mkdir -p tarfiles
    set tarfile = $homedir/tarfiles/$version.tar
endif

set PROJECTS       = `cat ${homedir}/PROJECTS`
set RUNS           = `cat ${homedir}/RUNS`
set DIRS           = `cat ${homedir}/DIRS`

# Initialize the number of saved and missing files and/or directories

set n_missing_dirs  = 0
set n_missing_files = 0

set missing_files   = ()
set missing_dirs    = ()

set n_saved_projects = 0
set saved_projects   = ()

set n_saved_runs     = 0
set saved_runs       = ()

set n_saved_dirs     = 0
set saved_dirs       = ()

set n_saved_files    = 0
set saved_files      = ()

set n_saved_bmk_dirs = 0
set saved_bmk_dirs   = ()

set files = `cat .files`

echo "========================================================================"
echo "Creating ${tarfile}"
echo "========================================================================"
echo "Source code        : *.py, *.c"
echo "------------------------------------------------------------------------"
echo "Files:             : $files"
echo "------------------------------------------------------------------------"
echo "papers             : `ls papers`"
echo "------------------------------------------------------------------------"
echo "C-shell scripts    : ./scripts/*.s"
echo "------------------------------------------------------------------------"
echo "runs               : $RUNS"
echo "projects           : $PROJECTS"
echo "========================================================================"

echo "Adding source code and shared object file"
echo "------------------------------------------------------------------------"

tar cvf $tarfile *.py *.c *.so .author --warning=no-file-changed >& /dev/null

echo "Adding files"
echo "------------------------------------------------------------------------"

tar rf $tarfile $files --warning=no-file-changed >& /dev/null

echo "Adding papers"

tar rf $tarfile papers --warning=no-file-changed

echo "Adding scripts"

tar rf $tarfile scripts --warning=no-file-changed

# Loop over projects

foreach project ($PROJECTS)

    echo "========================================================================"
    set project_dir = ./projects/$project/
    if (!(-e $project_dir)) then
        echo "Project directory does not exist, skipping: $project"
        @ n_missing_dirs = $n_missing_dirs + 1
        set missing_dirs =  ($missing_dirs $project_dir)
    else
        echo "Adding project: $project"
        @ n_saved_projects = $n_saved_projects + 1
        set saved_projects = ($saved_projects $project)

#       Loop over dirs to check for bmk directory

        foreach dir ($DIRS)
            if($dir == "bmk") then
                echo "------------------------------------------------------------------------"
                set dirname = ${project_dir}${dir}
                if(-e $dirname) then
                    echo "    Adding directory: $dirname"
                    @ n_saved_dirs = $n_saved_dirs + 1
                    set saved_drs =  ($saved_dirs $dirname)
                    tar rf $tarfile $dirname --warning=no-file-changed
                else
                    echo "    Directory does not exist, skipping: $dirname"
                    @ n_missing_dirs = $n_missing_dirs + 1
                    set missing_dirs = ($missing_dirs $dirname)
                endif
            endif
        end

#       Loop over runs per the RUNS file

        foreach run ($RUNS)

#           Each line in RUNS is "PROJECT RUN1 ... RUN" so ignore run = project (ie. first element of a line line in RUNS)

            if($run != $project) then
                echo "------------------------------------------------------------------------"
                set run_dir = ./projects/$project/$run/
                set in_file = ${project_dir}$run.in
                if (!(-e ${in_file})) then
                    echo "Input file not exist, skipping run: $run"
                    @ n_missing_files = $n_missing_files + 1
                    set missing_files = ($missing_files $in_file)
                else
                    echo "    Adding run: $run"
                    @ n_saved_runs = $n_saved_runs + 1
                    set saved_runs = ($saved_runs $run)
                    foreach filetype (in out bmk)
                        set filename = ${project_dir}$run.$filetype
                        if (-e $filename) then
                            echo "        Adding file     : $filename"
                            @ n_saved_files = $n_saved_files + 1
                            set saved_files = ($saved_files $in_file)
                            tar rf $tarfile  $filename --warning=no-file-changed
                        else
                            echo "        File does not exist, skipping: $filename"
                            @ n_missing_files = $n_missing_files + 1
                            set missing_files = ($missing_files $filename)
                        endif
                    end
                    foreach dir ($DIRS)
                        if($dir != "bmk") then
                            set dirname = ${project_dir}$run/$dir
                            if(-e $dirname) then
                                echo "        Adding directory: $dirname"
                                @ n_saved_dirs = $n_saved_dirs + 1
                                set saved_drs =  ($saved_dirs $dirname)
                                tar rf $tarfile $dirname --warning=no-file-changed
                            else
                                echo "        Directory does not exist, skipping: $dirname"
                                @ n_missing_dirs = $n_missing_dirs + 1
                                set missing_dirs = ($missing_dirs $dirname)
                            endif
                        endif
                    end
                endif
            endif
        end
echo "========================================================================"
echo "Saved $n_saved_runs run(s) from $n_saved_projects project(s)"
echo "------------------------------------------------------------------------"
echo "Saved a total of $n_saved_files files and $n_saved_dirs directories"
echo "========================================================================"
if($n_missing_files > 0 || $n_missing_dirs > 0) then
    echo "WARNING: Missing ${n_missing_files} files and ${n_missing_dirs} directories"
    echo "         These items have not been saved!"
    echo "========================================================================"
    echo "To resolve this, check that the PROJECTS, RUNS and DIRS files only specify"
    echo "the projects, runs and directories that exist and which you wish to save"
    echo "------------------------------------------------------------------------"
    echo "If .in files for a run exist but the corresponding .out & .bmk files plus the"
    echo "plots & data directories are missing, it probably means the run wasn't done"
    echo "------------------------------------------------------------------------"
    echo "If only *_bmk/ directories, it is not a problem as such directories are"
    echo "only create directories after running plot_bmk.py to generate benchmark"
    echo "plots such as scaling of Nu and Vrms vs Ra and corrsponding time series"
    echo "========================================================================"
    echo "Missing files"
    echo "------------------------------------------------------------------------"
    echo ${missing_files}
    echo "========================================================================"
    echo "Missing directories"
    echo "------------------------------------------------------------------------"
    echo ${missing_dirs}
    echo "========================================================================"
endif
echo "Tarfile created: "`du -s -h ${tarfile}`
echo "========================================================================"

popd >& /dev/null
