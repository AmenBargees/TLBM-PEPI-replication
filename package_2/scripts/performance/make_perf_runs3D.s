#!/usr/bin/csh

#==========================================================================
#
# Make the input and script files for 3D runs
# Note that the "" are NECESSARY where they are indicated
#
# Usage:
#
# csh make_perf_runs3D.s arg1 arg2 arg3 arg4 arg5 arg6 arg7 arg8
#
# csh perf_runs.s computer_name computer_type corespernode "dimensions" "scalings" "sizes" "nps" option
#
# where
#
# arg1 = computer_name = computer or cluster name
# arg2 = computer_type = cluster or workstation
#      (workstation runs are sequential, cluster scripts all be sent at once)
# arg3 = corespernode  = cores per node for cluster runs
# arg4 = dimensions    = must be "3D" as this script is for 3D runs
# arg5 = scalings      = "weak strong" for both weak & strong scaling runs or just "weak" or "strong"
# arg6 = sizes         = size of each axis relative to the axis lengths
# arg7 = nps           = core list such as "1 2 ... 512"   (for np    = 1, 2, ..., 512)
# arg8 = option        = option on how to modify dimensions for weak scaling
#                        where 1 = retain aspect ratio, 2 = (1 or 1.5) x log2 axes
#==========================================================================

set computer     = $1
set what         = $2
set corespernode = $3
set Ds           = "$4"
set scalings     = "$5"
set sizes        = "$6"
set nps          = "$7"
set option       = $8

set inputdir = ./perf_input/
set runsdir = ./runs_perf/
set outdir  = ./out_perf/
set tmpdir  = ./tmp_perf/
if($what == cluster) then
set jobsdir  = ./jobs_perf/
endif

#==========================================================================

# Set the base domain size for weak scaling runs plus time steps and inc

set nx0 = 32
set ny0 = 16
set nz0 = 16
# set nt0 = 50000 ORIGINAL ooo
# set inc0 = 5000
set nt0 = 2000
set inc0 = 200

# Set the base domain size for strong scaling runs plus time steps and inc

set nx1 = 64
set ny1 = 64
set nz1 = 64
set nt1 = 1000
set inc1 = 100

# Loop over base input files to add the correct options and nx,nz,nt for each run 
#
# Note that 30.txt is a base 3D input file in $inputdir
#
# 3D weak   scaling run performance files will start with 30 and have a zero separating the options
# 3D strong scaling run performance files will start with 31 and have a zero separating the options
# eg. 30010128 means a weak   scaling run with nx, ny and nz at size = 1 times the base domain and 128 cores
# eg. 31020256 means a strong scaling run with nx, ny and nz at size = 2 times the base domain and 256 cores
# The output performance files will be in $tmpdir//performance_${computer}_${scaling}_${D}_${size}_${np}

foreach D ( $Ds )

if($D == 3D) then

foreach scaling ($scalings)

foreach file0 (30)

# Loop over the domain size factor relative to base domain size of nx0 x ny0 x nz0
# Namely: (nx x ny x nz) = (nx0 x size) x (ny0 x size) x (nz0 x size) where size = 1, 2 etc
# nt is set such that runs take about the same time of a few minutes

foreach size ($sizes)

# Loop over how many cores (note: You need to be sure there
# is an appropriate if block for the cores you run with)
# Also, if you modify to add a new core count, you also need
# to modify the foreach np loop in the strong scaling runs below

foreach np ($nps)

#==============================================================================================
# Create a weak scaling input file and slurm script file 
#==============================================================================================

if($scaling == weak) then

@ file  = $file0

cp ${inputdir}/${file0}.txt ${runsdir}/${file}0${size}0${np}.txt

echo "performance_file = ${computer}_${scaling}_${D}_${size}_${np}" >> ${runsdir}/${file}0${size}0${np}.txt
echo "performance_append = False" >> ${runsdir}/${file}0${size}0${np}.txt
echo "tmp_dir = ${tmpdir}"  >> ${runsdir}/${file}0${size}0${np}.txt

# Get dimensions for this size and np

set nxnynz = "`${inputdir}/make_domain_size $D $np $size $nx0 $ny0 $nz0 $option`"
# last argument = 1 ... axis sizes that retain the aspect ratio
# last argument = 2 ... log2 axis sizes, or 1.5 x log2 axis sizes

set ix = 0
foreach n ($nxnynz)
if($ix == 0) then
set nx = $n
else if($ix == 1) then
set ny = $n
else if($ix == 2) then
set nz = $n
endif
@ ix = $ix + 1
end

@ nt  = $nt0 /  ($size * $size * $size)
@ inc = $inc0 / ($size * $size * $size)

if($nt < 1) then
@ nt = 1
endif
if($inc < 1) then
@ inc = 1
endif

echo "nx = "${nx} >> ${runsdir}/${file}0${size}0${np}.txt
echo "ny = "${ny} >> ${runsdir}/${file}0${size}0${np}.txt
echo "nz = "${nz} >> ${runsdir}/${file}0${size}0${np}.txt
echo "nt = "${nt} >> ${runsdir}/${file}0${size}0${np}.txt
echo "Inc_timing  = "$inc >> ${runsdir}/${file}0${size}0${np}.txt
echo "Inc_print_U = "$inc >> ${runsdir}/${file}0${size}0${np}.txt

echo "=============================================="
echo "Weak scaling 3D input file created: "${runsdir}${file}0${size}0${np}.txt
echo "size = "$size
echo "np = "${np}
echo "nx = "${nx}
echo "ny = "${ny}
echo "nz = "${nz}
echo "nt = "${nt}
@ n = ${nx} * ${ny} * ${nz} / ${np}
echo "n/core = " $n
@ dt = ($nt * $n) / 100000
echo "Relative calculation time = " $dt

#===================================================
# Change the slurm sbatch commands as necessary
#===================================================

if($what == cluster) then

set batchscript = ${runsdir}/run_${computer}_${scaling}_${D}_${size}_${np}.s

echo "----------------------------------------------"
echo "Creating batch file for the run file:"
echo "${runsdir}${file}0${size}0${np}.txt"

@ nodes = $np / $corespernode
@ nrem = $np % $corespernode
if($nrem > 0) then
@ nodes = $nodes + 1
endif
echo "#\!/usr/bin/csh" > ${batchscript}
echo "#SBATCH --job-name=LBM" >> ${batchscript}
echo "#SBATCH --output=${jobsdir}/job%j.out" >> ${batchscript}
echo "#SBATCH --error=${jobsdir}/job%j.err" >> ${batchscript}
echo "#SBATCH --nodes="$nodes >> ${batchscript}
echo "#SBATCH --partition=DRP" >> ${batchscript}
echo "#SBATCH --ntasks-per-node=1" >> ${batchscript}
echo "#SBATCH --cpus-per-task="${corespernode} >> ${batchscript}
echo "#SBATCH --exclusive" >> ${batchscript}
echo "#SBATCH --mem=0" >> ${batchscript}
echo "module load intel/Python-3.6.4" >> ${batchscript}
echo "source ~/ctypes/ctypes-1.0.2" >> ${batchscript}

echo "mpiexec -np ${np} python3 TLBM_fast.py ${file}0${size}0${np} ${runsdir} > & ${outdir}/${file}0${size}0${np}.out" >> ${batchscript}

echo "Batch file created: "${batchscript}
echo "----------------------------------------------"

endif

endif

#==============================================================================================
# Create a strong scaling input file and slurm script file 
#==============================================================================================

if($scaling == strong) then

@ file  = $file0 + 1

cp ${inputdir}/${file0}.txt ${runsdir}/${file}0${size}0${np}.txt

echo "performance_file = ${computer}_${scaling}_${D}_${size}_${np}" >> ${runsdir}/${file}0${size}0${np}.txt
echo "performance_append = False" >> ${runsdir}/${file}0${size}0${np}.txt
echo "tmp_dir = ${tmpdir}"  >> ${runsdir}/${file}0${size}0${np}.txt

@ nx = $nx1 * $size
@ ny = $ny1 * $size
@ nz = $nz1 * $size
@ nt = $nt1 * $np / ($size * $size * $size)
@ inc = $inc1 * $np / ($size * $size * $size)

if($nt < 1) then
@ nt = 1
endif
if($inc < 1) then
@ inc = 1
endif

echo "nx = "$nx >> ${runsdir}/${file}0${size}0${np}.txt
echo "ny = "$ny >> ${runsdir}/${file}0${size}0${np}.txt
echo "nz = "$nz >> ${runsdir}/${file}0${size}0${np}.txt
echo "nt = "$nt >> ${runsdir}/${file}0${size}0${np}.txt
echo "Inc_timing  = "$inc >> ${runsdir}/${file}0${size}0${np}.txt
echo "Inc_print_U = "$inc >> ${runsdir}/${file}0${size}0${np}.txt

echo "=============================================="
echo "Strong scaling 3D input file created: "${runsdir}${file}0${size}0${np}.txt
echo "size = "$size
echo "np = "${np}
echo "np = "${np}
echo "nx = "${nx}
echo "nz = "${nz}
echo "nt = "${nt}
@ n = ${nx} * ${ny} * ${nz} / ${np}
echo "n/core = " $n
@ dt = ($nt * $n) / 100000
echo "Relative calculation time = " $dt

#===================================================
# Change the slurm sbatch commands as necessary
#===================================================

if($what == cluster) then

set batchscript = ${runsdir}/run_${computer}_${scaling}_${D}_${size}_${np}.s

echo "----------------------------------------------"
echo "Creating batch file"

@ nodes = $np / $corespernode
@ nrem = $np % $corespernode
if($nrem > 0) then
@ nodes = $nodes + 1
endif
echo "#\!/usr/bin/csh" > ${batchscript}
echo "#SBATCH --job-name=LBM" >> ${batchscript}
echo "#SBATCH --output=${jobsdir}/job%j.out" >> ${batchscript}
echo "#SBATCH --error=${jobsdir}/job%j.err" >> ${batchscript}
echo "#SBATCH --nodes="$nodes >> ${batchscript}
echo "#SBATCH --partition=DRP" >> ${batchscript}
echo "#SBATCH --ntasks-per-node=1" >> ${batchscript}
echo "#SBATCH --cpus-per-task="${corespernode} >> ${batchscript}
echo "#SBATCH --exclusive" >> ${batchscript}
echo "#SBATCH --mem=0" >> ${batchscript}
echo "module load intel/Python-3.6.4" >> ${batchscript}
echo "source ~/ctypes/ctypes-1.0.2" >> ${batchscript}

echo "mpiexec -np ${np} python3 TLBM_fast.py ${file}0${size}0${np} ${runsdir} > & ${outdir}/${file}0${size}0${np}.out" >> ${batchscript}

echo "Batch file created: "${batchscript}
echo "----------------------------------------------"

endif

endif

end
end

end
end

endif

end
