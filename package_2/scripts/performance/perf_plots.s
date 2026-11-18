#!/bin/csh -f

#-------------------------------------------------------------------------------------------------
# Generate scaling and gflops plots in directory ./results_perf
#
# First run perf_runs.s to generate timing files in ./tmp_perf
# Then run perf_plots.s as described below. Note that the ""
# are NECESSARY where they are indicated
#
# Usage:
#
# csh perf_plots.s computer_name dimensions scalings sizes nps
#
# or
#
# csh perf_plots.s arg1 arg2 arg3 arg4 arg5
#
# arg1 = computer_name = computer or cluster name
# arg2 = dimensions    = "2D 3D" for both 2D and 3D scaling runs, or just "2D" or "3D"
# arg3 = scalings      = "weak strong" for both weak & strong scaling runs or just "weak" or "strong"
# arg4 = sizes         = size of each axis relative to the axis lengths
# arg5 = nps           = core list such as "1 2 ... 512"   (for np    = 1, 2, ..., 512)
#-------------------------------------------------------------------------------------------------

set computer     = $1
set Ds           = "$2"
set scalings     = "$3"
set sizes        = "$4"
set nps          = "$5"


# Set directories for input files, run files, output files and slurm output files

set inputdir = ./perf_input/
set runsdir = ./runs_perf/
set outdir  = ./out_perf/
set tmpdir  = ./tmp_perf/
set resultsdir  = ./results_perf/
mkdir -p ${resultsdir}

foreach D ( $Ds )
foreach scaling ( $scalings )
foreach size ( $sizes )

rm -f ${resultsdir}/performance_${computer}_${scaling}_${D}_${size}.txt

foreach np ( $nps )

cat ${tmpdir}/performance_${computer}_${scaling}_${D}_${size}_${np}.txt >> ${resultsdir}/performance_${computer}_${scaling}_${D}_${size}.txt

end

echo "Concatenated performance files to generate "performance_${computer}_${scaling}_${D}_${size}.txt

end
end
end

foreach D ( $Ds )
foreach scaling ( $scalings )
foreach size ( $sizes )

python3 performance.py $computer $D $scaling $size

end
end
end
