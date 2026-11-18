
# Set graphics backend to TkAgg so Python can create and save figures without an X server

setenv MPLBACKEND TkAgg

# Useful abbreviations for common commands

alias  H history
alias  P pushd
alias  p popd
alias  M more
alias  s source
alias  j jobs
alias  k  kill -9

# Set history to be long and unlimit limits

set    history=1000000
unlimit

# Set GTXLBM home directory and script directory, add script directory to PATH and rehash

setenv homedir $HOME/src/tlbm
alias  home pushd $homedir

setenv scriptdir $homedir/scripts
setenv PATH ./:${scriptdir}:$PATH

rehash

# Alias for developers to update the GTXLBM version and to compile GTXLBM.c

alias  update $scriptdir/update.s

# Alias to get IP address and date
 
alias  ip_address ip_address.s
alias  get_date get_date.s

# Aliases to do GTXLBM run(s) for specified run(s) and project(s)

alias  run     $scriptdir/run.s
alias  runs    $scriptdir/runs.s
alias  run_ex  $scriptdir/run_ex.s
alias  runb    $scriptdir/runb.s

# Abbreviations to ps and grep for mpiexec, python3, GTXLBM, and for GTXLBM and plot_npy

alias  psmpi  ps -alF \| grep mpiexec
alias  psm    ps -alF \| grep mpiexec
alias  pspy   ps -alF \| grep python3
alias  psmr   ps -alF \| grep mpiexec \| grep GTXLBM_run.py \| awk \'\{print \$4\}\'
alias  pspr   ps -alF \| grep python3 \| grep GTXLBM_run.py \| awk \'\{print \$4\}\'
alias  psr    ps -alF \| grep mpiexec \| grep GTXLBM_run.py \| grep \`cat $homedir/PROJECT\` \| grep \`cat $homedir/RUN\` \| awk \'\{print \$4\}\'
alias  psp    ps -alF \| grep python3 \| grep plot_npy
alias  pst    ps -alF \| grep GTXLBM

# Abbreviations to kill mpiexec, python3, and GTXLBM runs 

alias  kmpi k \`psmpi\`
alias  kpy  k \`pspy\`
alias  kmr  k \`psmr\`
alias  kpr  k \`pspr\`
alias  kr   k \`psr\`
alias  kp   k \`psp\`

# Aliases to generate tar files and backup directories of the GTXLBM home directory and run(s) and project(s)

alias  t      source $scriptdir/tar.s     # Abbreviation for tar.s
alias  tarp   tarp.s                      # Abbreviation for tarp.s
alias  tarall tarall.s                    # Abbreviation for tarall.s
alias  bak    source $scriptdir/bak.s     # Abbreviation for bak.s
alias  b      source $scriptdir/bak.s     # Abbreviation for bak.s
alias  bakp   bakp.s                      # Abbreviation for bakp.s
alias  bakall bakall.s                    # Abbreviation for bakall.s

# Alias to display (png) files using eog without warning messages

alias  eg $homedir/scripts/eg.s           # Abbreviation for eog

# Aliases to display plot (png) files from GTXLBM run of a project

alias  ep $homedir/scripts/ep.s           # Abbreviation for eog      (best choice as no X server is needed so it will always work)
alias  dp $homedir/scripts/dp.s           # Abbreviation for display  (good choice if you have an X server)
alias  Cp $homedir/scripts/Cp.s           # Abbreviation for cixel    (not recommended as it overwrites the terminal & the image is low resoluton)
alias  Ip $homedir/scripts/Ip.s           # Abbreviation for im2sixel
alias  fp $homedir/scripts/fp.s           # Abbreviation for fbi

# Aliases to movie players

alias  mplayer  mplayer -vo x11           # Abbreviation for mplayer (good option, requires an X server)
alias  mpl      mplayer                   # Abbreviation for mplayer
alias  vplayer  vlc                       # Abbreviation for vlc     (highly recommended as the movie window can be resized and paused etc)
alias  vpl      vlc                       # Abbreviation for vlc

# Aliases to pushd to GTXLBM directories

alias  PP source $scriptdir/PP.s          # Abbreviation to pushd to the project directory
alias  pb source $scriptdir/pb.s          # Abbreviation to pushd to the bmk     directory of a project
alias  pp source $scriptdir/pp.s          # Abbreviation to pushd to the plots   directory of a run of a project
alias  pd source $scriptdir/pd.s          # Abbreviation to pushd to the data    directory of a run of a project
alias  pm source $scriptdir/pm.s          # Abbreviation to pushd to the movies  directory of a run of a project

# Aliases to ls plot, data and movie directories, or files within these directories

alias  lsP $homedir/scripts/lsP.s         # Abbreviation to ls the project directory
alias  lsb $homedir/scripts/lsb.s         # Abbreviation to ls the bmk     directory of of a project
alias  lsp $homedir/scripts/lsp.s         # Abbreviation to ls the plots   directory of a run of a project
alias  lsd $homedir/scripts/lsd.s         # Abbreviation to ls the data    directory of a run of a project
alias  lsm $homedir/scripts/lsm.s         # Abbreviation to ls the movies  directory of a run of a project

# Aliases to vi input, output, error and bmk files of a GTXLBM_run.py run, plot generation
# output file of a plot_npy run, and the RUN, RUNS, PROJECT and PROJECTS files

alias  vit    vi $homedir/GTXLBM.py
alias  vip    vi $homedir/plot_npy.py
alias  vib    vi $homedir/plot_bmk.py
alias  vii    vi $homedir/projects/\`cat $homedir/PROJECT\`/\`cat $homedir/RUN\`.in
alias  vio    vi $homedir/projects/\`cat $homedir/PROJECT\`/\`cat $homedir/RUN\`.out
alias  vie    vi $homedir/projects/\`cat $homedir/PROJECT\`/\`cat $homedir/RUN\`.err
alias  vibmk  vi $homedir/projects/\`cat $homedir/PROJECT\`/\`cat $homedir/RUN\`.bmk
alias  vir    vi $homedir/RUN
alias  virs   vi $homedir/RUNS
alias  vip    vi $homedir/PROJECT
alias  vips   vi $homedir/PROJECTS

# Aliases to tail input, output, error and bmk files of a GTXLBM_run.py run, and plot generation output file of a plot_npy run

alias  ti  $scriptdir/ti.s                 # Abbreviation to tail the .out file of a GTXLBM run            of a project
alias  to  $scriptdir/to.s                 # Abbreviation to tail the .out file of a GTXLBM run            of a project
alias  te  $scriptdir/te.s                 # Abbreviation to tail the .err file of a GTXLBM run            of a project
alias  tb  $scriptdir/tb.s                 # Abbreviation to tail the .bmk file of a GTXLBM run            of a project
alias  tp  $scriptdir/tp.s                 # Abbreviation to tail the .out file of plot_npy run for run(s) of a project
alias  tr  tail $homedir/RUN
alias  trs tail $homedir/RUNS
alias  tP  tail $homedir/PROJECT
alias  tPs tail $homedir/PROJECTS

# Aliases to more input, output, error and bmk files of a GTXLBM_run.py run, and plot generation output file of a plot_npy run

alias  mi $scriptdir/mi.s                 # Abbreviation to more the .in  file of a GTXLBM run            of a project
alias  mo $scriptdir/mo.s                 # Abbreviation to more the .in  file of a GTXLBM run            of a project
alias  me $scriptdir/me.s                 # Abbreviation to tail the .err file of a GTXLBM run            of a project
alias  mb $scriptdir/mb.s                 # Abbreviation to tail the .err file of a GTXLBM run            of a project
alias  mp $scriptdir/mp.s                 # Abbreviation to tail the .out file of plot_npy run for run(s) of a project
alias  mr  tail $homedir/RUN
alias  mrs tail $homedir/RUNS
alias  mP  tail $homedir/PROJECT
alias  mPs tail $homedir/PROJECTS

# Aliases to grep input, output, error and bmk files of a GTXLBM_run.py run, and plot generation output file of a plot_npy run

alias  gi $scriptdir/gi.s                 # Abbreviation to grep the .in  file of a GTXLBM   run of a project
alias  go $scriptdir/go.s                 # Abbreviation to grep the .out file of a GTXLBM   run of a project
alias  ge $scriptdir/ge.s                 # Abbreviation to grep the .err file of a GTXLBM   run of a project
alias  gb $scriptdir/gb.s                 # Abbreviation to grep the .bmk file of a GTXLBM   run of a project
alias  gp $scriptdir/gp.s                 # Abbreviation to grep the .out file of a plot_npy run of a project

# Aliases to generate plot (png) files from data files (npy) of GTXLBM runs
 
alias  plot    source $scriptdir/plot.s
alias  plots   source $scriptdir/plots.s
alias  plot_ex source $scriptdir/plot_ex.s

# Aliases to generate movie files from plot (png) files generated by GTXLBM or plot_npy runs

alias  movie     source $scriptdir/movie.s
alias  movies    source $scriptdir/movies.s
alias  movies_ex source $scriptdir/movies_ex.s
alias  movie2D   source $scriptdir/movie2D.s
alias  movie3D   source $scriptdir/movie3D.s

# Alias to clean files and directories of GTXLBM run by removing the .out, .bmk and .err files, and the plot, data and movie directories of the run of a project

alias  cleanr source $scriptdir/cleanr.s
