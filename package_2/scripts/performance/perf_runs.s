#!/usr/bin/csh -f

#---------------------------------------------------------
# Prepare input files and perform scaling runs
#
# Usage:
#
# csh perf_runs.s computer_name computer_type corespernode "dimensions" "scalings" "sizes" "nps" option
#
# ie.
#
# csh perf_runs.s arg1 arg2 arg3 arg4 arg5 arg6 arg7 arg8
#
#
# where
#
# arg1 = computer_name = computer or cluster name
# arg2 = computer_type = cluster or workstation
#      (workstation runs are sequential, cluster scripts all be sent at once)
# arg3 = corespernode  = cores per node for cluster runs
# arg4 = dimensions    = "2D 3D" for both 2D and 3D scaling runs, or just "2D" or "3D"
# arg5 = scalings      = "weak strong" for both weak & strong scaling runs or just "weak" or "strong"
# arg6 = sizes         = size of each axis relative to the axis lengths
# arg7 = nps           = core list such as "1 2 ... 512"   (for np    = 1, 2, ..., 512)
# arg8 = option        = option on how to modify dimensions for weak scaling
#                        where 1 = retain aspect ratio, 2 = (1 or 1.5) x log2 axes
#---------------------------------------------------------

set computer     = $1
set what         = $2
set corespernode = $3
set Ds           = "$4"
set scalings     = "$5"
set sizes        = "$6"
set nps          = "$7"
set option       = "$8"

# Set directories for input files, run files, output files and slurm output files

set inputdir = ./perf_input/
set runsdir  = ./runs_perf/
set outdir   = ./out_perf/
set tmpdir   = ./tmp_perf/

# Make directories if needed

mkdir -p $runsdir
mkdir -p $outdir
mkdir -p $tmpdir

if($what == cluster) then
set jobsdir  = ./jobs_perf/
mkdir -p $jobsdir
endif

# Compile make_domain_size

gcc ${inputdir}/make_domain_size.c -lm -o ${inputdir}/make_domain_size

echo ""
echo "==========================================="
echo "Making input files and scripts as necessary"
echo "==========================================="
echo ""
sleep 1

foreach D ( $Ds )
if($D == 2D) then
echo ""
echo "====================="
echo "Making 2D perf files"
echo "====================="
echo ""
csh ${inputdir}/make_perf_runs2D.s $computer $what $corespernode "$D" "$scalings" "$sizes" "$nps" $option
else if($D == 3D) then
echo ""
echo "====================="
echo "Making 3D perf files"
echo "====================="
echo ""
csh ${inputdir}/make_perf_runs3D.s $computer $what $corespernode "$D" "$scalings" "$sizes" "$nps" $option
endif
end

echo ""

if($what == cluster) then

echo "================================================"
echo "All input files and slurm scripts have been made"
echo "================================================"

else if($what == workstation) then

echo "=============================="
echo "All input files have been made"
echo "=============================="

endif

echo ""
sleep 1

echo ""
if($what == cluster) then

echo "============================="
echo "Submitting all sbatch scripts"
echo "============================="

else if($what == workstation) then

echo "============================"
echo "Starting runs on workstation"
date
echo "============================"

endif

echo ""
sleep 1


# Loop for 2D and 3D tests

foreach D ( $Ds ) 

# Loop for weak and/or strong scaling runs

foreach scaling ( $scalings )

# Loop for sizes

foreach size ( $sizes)

# Loop for number of cores

foreach np ( $nps )

if ($what == cluster) then

echo "Starting batch run: ./${runsdir}/"run_${computer}_${scaling}_${D}_${size}_${np}.s

sbatch --exclusive ./${runsdir}/run_${computer}_${scaling}_${D}_${size}_${np}.s

else if ($what == workstation) then

if($D == 2D) then
@ file = 20
else if($D == 3D) then
@ file = 30
endif
if($scaling == strong) then
@ file = $file + 1
endif

echo "====================================================================="
echo "Doing run "${file}0${size}0${np}
echo "D       = "$D
echo "Scaling = "$scaling
echo "Size    = "$size
echo "np      = "$np
echo "---------------------------------------------------------------------"
echo "Started: mpiexec -np ${np} python3 TLBM.py ${file}0${size}0${np} ${runsdir} > & ${outdir}${file}0${size}0${np}.out"
mpiexec -np ${np} python3 TLBM.py ${file}0${size}0${np} ${runsdir} > & ${outdir}${file}0${size}0${np}.out
echo "Finished: mpiexec -np ${np} python3 TLBM.py ${file}0${size}0${np} ${runsdir} > & ${outdir}${file}0${size}0${np}.out"

endif

end
end
end
end

if($what == cluster) then

echo ""
echo "============================"
echo "All sbatch scripts submitted"
echo "============================"
echo ""

else if($what == workstation) then

echo "====================="
echo "ALL RUNS ARE COMPLETE"
date
echo "====================="

endif

