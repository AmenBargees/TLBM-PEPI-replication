#!/bin/csh -f


set PROJECTS = `cat $homedir/PROJECTS`
set RUNS     = `cat $homedir/RUNS`

if($#argv == 1 && $1 == "-h") then
    echo "------------------------------------------------------------------------"
    echo "Foreground run(s) of GTXLBM for multiple run(s) & project(s)"
    echo "------------------------------------------------------------------------"
    echo "Usage   : runs.s [cores]"
    echo "------------------------------------------------------------------------"
    echo "Examples: runs"
    echo "          runs 192"
    echo "------------------------------------------------------------------------"
    echo "Notes: (1) Run(s) & project(s) per RUN and PROJECT files"
    echo "       (2) Cores per CORES file if unspecified"
    echo "       (3) Type tr or tr.s to tail the .out file & check on run status"
    echo "       (4) Requires the run input file: ./projects/project/run.in"
    echo "       (5) Output files & directories :"
    echo " "
    echo "           ./projects/project/run.out"
    echo "           ./projects/project/run/data/*.npy"
    echo "           ./projects/project/run/plots/*.png"
    echo "           ./projects/project/run/movies/"
    echo "------------------------------------------------------------------------"
    exit
endif

if($#argv == 1) then
    set cores = $1
else
    set cores = `cat $homedir/CORES`
endif

echo "=============================================================================="
echo "Starting foreground run(s) of GTXLBM on ${cores} cores."
echo "To check run status, type sr or sr.s"
echo "------------------------------------------------------------------------------"
echo "runs     : $RUNS"
echo "projects : $PROJECTS"

pushd $homedir >& /dev/null

set n_missing_projects = 0
set missing_projects   = ()

set n_projects         = 0
set projects           = ()

set n_missing_runs     = 0
set missing_runs       = ()
set missing_in_files   = ()

set n_runs             = 0
set run_out_files      = ()

foreach project ($PROJECTS)
    set i_project      = 0
    set n_runs_project = 0
    set end_project        = false
    set project_dir = $homedir/projects/$project
    echo "=============================================================================="
    echo "Doing runs for project $project"
    echo "=============================================================================="
    if(!(-e $project_dir)) then
        echo "Project directory does not exist: $project_dir"

        @ n_missing_projects = $n_missing_projects + 1
        set missing_projects = ($missing_projects  $project_dir)

        @ n_projects = $n_projects + 1
        set projects = ($projects  $project)
    else
        foreach run ($RUNS)
            if($run == $project) then
                set n_runs_project = 0
                set end_project    = false
            else
                set next_project_dir = $homedir/projects/$run
                set run_in_file      = $homedir/projects/$project/$run.in
                if(-e $next_project_dir) then
                    set end_project = true
                else
                    if(!(-e $run_in_file)) then
                        echo "Run could not be done for non-existent input file:"
                        echo "    $run_in_file"
                        @ n_missing_runs     = $n_missing_runs + 1
                        set missing_runs     = ($missing_runs     $run    )
                        set missing_in_files = ($missing_in_files $run_in_file)
                    else if (-e $run_in_file) then
                        @ n_runs_project = $n_runs_project + 1
                        echo "Preparing to start run for input file: ${run_in_file}"
                        echo "    File & directory cleanup ..."
                        echo -n "        rm $homedir/projects/$project/$run.out ... "
                        rm $homedir/projects/$project/$run.out        >& /dev/null
                        echo "Done"
        
                        echo -n "        rm $homedir/projects/$project/$run.bmk ... " 
                        rm $homedir/projects/$project/$run.bmk        >& /dev/null
                        echo "Done"

                        echo -n "        rm -r $homedir/projects/$project/$run/plots/ ... "
                        rm -r $homedir/projects/$project/$run/plots/  >& /dev/null
                        echo "Done"
    
                        echo -n "        rm -r $homedir/projects/$project/$run/data/ ... "
                        rm -r $homedir/projects/$project/$run/data/   >& /dev/null
                        echo "Done"

                        echo -n "        rm -r $homedir/projects/$project/$run/movies/ ... "
                        rm -r $homedir/projects/$project/$run/movies/ >& /dev/null
                        echo "Done"     
    
                        echo -n "Starting run on $cores cores: "
                        date
                        echo "Output file : $homedir/projects/$project/$run.out"
                        echo -n "mpiexec -np $cores python3 GTXLBM_run.py $project $run >& $homedir/projects/$project/$run.err ... running ... "
                        mpiexec -np $cores python3 GTXLBM_run.py $project $run >& $homedir/projects/$project/$run.err
                        echo "complete!"
                        @ n_runs = $n_runs + 1
                        set run_out_files = ($run_out_files $homedir/projects/$project/$run.out)
                    endif
                endif
            endif
        end
    endif
end

echo "=============================================================================="
echo " "
echo " "
echo "=============================================================================="
echo "All runs completed!"
echo "=============================================================================="
echo "SUMMARY"
echo "=============================================================================="
if($n_missing_projects > 0 || $n_missing_runs > 0) then
    echo "$n_missing_runs runs could not be done due to non-existent input files!"
    echo "------------------------------------------------------------------------------"
    if($n_missing_projects > 0) then
        if($n_missing_projects == 1) then
            echo "$n_missing_projects missing project directory:"
        else
            echo "$n_missing_projects missing project directories:"
        endif
        foreach i ($missing_projects)
            echo "    $i"
        end
        if($n_missing_runs > 0) then
            echo "------------------------------------------------------------------------------"
        endif
    endif
    if($n_missing_runs > 0) then
        if($n_missing_runs == 1) then
            echo "$n_missing_runs missing input file:"
        else
            echo "$n_missing_runs missing input files:"
        endif
        foreach i ($missing_in_files)
            echo "    $i"
        end
    endif
    if($n_runs > 0) then
        echo "=============================================================================="
    endif
endif

if($n_runs > 0) then
    if($n_runs == 1) then
        echo "$n_runs run completed! Output file:"
    else
        echo "$n_runs runs completed! Output files:"
    endif
    echo "------------------------------------------------------------------------------"
    foreach i ($run_out_files)
        echo "$i"
    end
endif
echo "=============================================================================="
popd >& /dev/null
