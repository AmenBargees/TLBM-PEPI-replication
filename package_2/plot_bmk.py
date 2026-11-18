#================================================================================
# Set suite versions, suites, BCs, and what time series
# and scaling plots to make, and set figsize below here
#================================================================================

#--------------------------------------------------------------------------
# Start of input options

# Set paper to PEPI-method or GJI-benchmark to get plot options for the figures in these papers

# paper         = "PEPI-method"
paper         = "GJI-benchmark"

project       = "benchmarks_annulus"
Ras_suite     = [4,5,6,7,8,9]
whats_suite   = ["Pr=1000"]

if paper == "PEPI-method":
    codes         = ["GTXLBM"]
    BCs           = ["free_slip"]
else:
    codes         = ["ASPECT","GTXLBM"]
    BCs           = ["free_slip","no_slip"]

test          = False ; resolution    = "" # Resolution allows allows different tests -> different png files. e.g. resolution = "A"

whats_tseries = ["Nu","Vrms","TBLW","Ma"]
whats_scaling = ["Nu","Vrms"]

averaging_window_start = 0.6 ; averaging_window_end   = 1

# End of options
#--------------------------------------------------------------------------

import numpy as np
import sys

nruns_max     = max(10,len(Ras_suite))
TBLW_data = True
fouts     = [sys.stderr,sys.stdout]
whats_version = [""]        # Always leave this empty!

Ras_ASPECT  = 10**np.array([4   ,5   ,6   ,7    ,8    ,9    ],dtype=float)
Nus_ASPECT_no_slip  = np.array([2.437,4.549,7.477,13.59,26.71,52.78],dtype=float)
Vrms_ASPECT_no_slip = np.array([19.3,79.4,289.1,1184,4323,18470],dtype=float)
Nus_ASPECT_free_slip  = np.array([4.62,9.41,18.7,39.1,72.5,149.8],dtype=float)
Vrms_ASPECT_free_slip = np.array([46.7,196,799,3355,13890,57270],dtype=float)

plot_GL_Iu              = False

# Plot fontsizes and option to plot GL Iu theory on scaling plots

figsize                    = "small"

f_line   = 1.8 # Factor wider lines and larger markers for ASPECT scaling plots
f_marker = 1.5 # Factor wider lines and larger markers for ASPECT scaling plots
if figsize == "small":
    left    = 0.11 ; bottom    = 0.12 ; top    = 0.93 ; right    = 0.975
    left_sc = 0.12 ; bottom_sc = 0.13 ; top_sc = 0.93 ; right_sc = 0.98
    linewidth               = 2.5
    if (paper == "PEPI-method" or paper == "GJI-benchmark") and not test:
        if paper == "PEPI-method":
            markersize = 10
            linewidth = 3.5
        else:
            markersize = 8
            linewidth = 2.5
        xsize                   = 7
        ysize                   = 5
    else                     :
        markersize = 6
        xsize                   = 5.5
        ysize                   = 4

    fontsize_title          = 16
    fontsize_axis_labels    = 16
    fontsize_axis_tics      = 10
    fontsize_legend_tseries = 12
    if paper != "PEPI-method":
        fontsize_legend_scaling = 12
    else:
        fontsize_legend_scaling = 14
else:
    left = 0.09 ; bottom = 0.1 ; top = 0.95 ; right = 0.985
    linewidth               = 3
    markersize              = 12
    xsize                   = 10
    ysize                   = 7
    fontsize_title          = 14
    fontsize_axis_labels    = 14
    fontsize_axis_tics      = 12
    fontsize_legend_tseries = 14
    fontsize_legend_scaling = 14

# Import libraries and get GTXLBM parameters

#================================================================================
# Set i_X for the bmk file format where X is the index of each variable in the bmk file
#================================================================================

i_t_p_int       = 0
i_timescale_t_p = 1
i_Ra            = 2
i_Pr            = 3
i_A             = 4
i_Ma            = 5
i_Tbar          = 6

def get_indices(words):

    global i_Nu, i_Nu_fit, i_Vrms, i_Re, i_TBLW, i_TBLW_fit, i_Nu_u, i_Nu_d

    i_Nu           = 7
    if len(words) >= 11:
        i_Vrms     = 8
        i_Re       = 9
        i_TBLW     = 10
        i_Nu_u     = 11
        i_Nu_d     = 12
    else:
        i_Nu_fit   = 8
        i_Vrms     = 9
        i_Re       = 10
        i_TBLW     = 11
        i_TBLW_fit = 12

#================================================================================
# Import Python libraries and GTXLBM libraries
#================================================================================

import matplotlib
matplotlib.use("Agg") # This is because Agg works without an x server, and avoids problems should the x server
import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import axes3d, Axes3D
import sys
import os
from os.path import exists

LaTeX = False 
if LaTeX: matplotlib.rcParams['text.usetex'] = True

import utilities
from utilities import *
from utilities import get_lattice_viscosity
from utilities import get_annular_lattice
from utilities import get_R_and_L
from utilities import get_time_value_text_int_timescale_from_time_spec
from utilities import timesteps_to_time_text
from utilities import timetext
from utilities import timevalue_to_secs
import get_parameters 
from get_parameters import get_parameters

#================================================================================
# Function to plot the scaling relationship
#================================================================================

def plot_scaling(title,filename,x,y,xx,yy,xlab,ylab):

    global code, what, log_xx_min, log_xx_max, log_yy_min, log_yy_max

    error     = (y_data - y_fit)/y_fit
    rms_error = np.sqrt(np.sum(error**2)/len(y_fit))

    log_xx_min  = min(log_xx_min,np.log10(np.min(xx)))
    log_xx_max  = max(log_xx_max,np.log10(np.max(xx)))

    log_yy_min  = min(log_yy_min,np.log10(np.min(yy)))
    log_yy_max  = max(log_yy_max,np.log10(np.max(yy)))

    log_yy_min  = min(log_yy_min,np.log10(np.min(y )))
    log_yy_max  = max(log_yy_max,np.log10(np.max(y )))

    if   BC == "no_slip"  : BClab = "zero-slip"
    elif BC == "free_slip": BClab = "free-slip"
    elif BC == "periodic" : BClab = "periodic"
    if code == "GTXLBM":
        if what == "Vrms":
            if paper == "PEPI-method" and len(BCs) == 1:
                labd = r'$\overline{V}_{rms}$ data'
                labf = r'$V_{rms} = \, %5.3f \, Ra^{%5.3f \pm %5.3f}$' %(const_fit,power_fit,stddev)
            else:
                labd = r'%s (%s): $\overline{V}_{rms}(Ra)$' %(code,BClab)
                labf = r'%s (%s): $V_{rms}(Ra) = \, %5.3f \, Ra^{%5.3f \pm %5.3f}$' %(code,BClab,const_fit,power_fit,stddev)
        else:
            if paper == "PEPI-method" and len(BCs) == 1:
                labd = r'$\overline{%s}$ data' %(what)
                labf = r'$%s = \, %5.3f \, Ra^{%5.3f \pm %5.3f}$' %(what,const_fit,power_fit,stddev)
            else:
                labd = r'%s (%s): $\overline{%s}(Ra)$' %(code,BClab,what)
                labf = r'%s (%s): $%s(Ra)      = \, %5.3f \, Ra^{%5.3f \pm %5.3f}$' %(code,BClab,what,const_fit,power_fit,stddev)

        if   BC == "no_slip":
            plt.plot( x, y,"o",color='blue',markersize=markersize,markerfacecolor='blue',markeredgecolor='black')
            plt.plot(xx,yy,"-",color='blue',label=labf,linewidth=linewidth)
        elif BC == "free_slip":
            if paper == "PEPI-method":
                plt.plot( x, y,"o",color='blue',label=labd,markersize=markersize,markerfacecolor='blue',markeredgecolor='black')
                plt.plot(xx,yy,"-",color='blue',label=labf,linewidth=linewidth)
            else:
                plt.plot( x, y,"o",color='red',markersize=markersize,markerfacecolor='red',markeredgecolor='black')
                plt.plot(xx,yy,"-",color='red',label=labf,linewidth=linewidth)

    elif code == "ASPECT":
        if what == "Vrms":
            if paper == "PEPI-method" and len(BCs) == 1:
                labd = r'%s  : $\overline{V}_{rms}(Ra)$' %(code)
                labf = r'%s  : $V_{rms}(Ra) = \, %5.3f \, Ra^{%5.3f}$' %(code,const_fit,power_fit)
            else:
                labd = r'%s  (%s): $\overline{V}_{rms}(Ra)$' %(code,BClab)
                labf = r'%s  (%s): $V_{rms}(Ra) = \, %5.3f \, Ra^{%5.3f \pm %5.3f}$' %(code,BClab,const_fit,power_fit,stddev)
        else:
            if paper == "PEPI-method" and len(BCs) == 1:
                labd = r'%s  : $\overline{%s}(Ra)$' %(code,what)
                labf = r'%s  : $%s(Ra)      = \, %5.3f \, Ra^{%5.3f}$' %(code,what,const_fit,power_fit)
            else:
                labd = r'%s  (%s): $\overline{%s}(Ra)$' %(code,BClab,what)
                labf = r'%s  (%s): $%s(Ra)      = \, %5.3f \, Ra^{%5.3f \pm %5.3f}$' %(code,BClab,what,const_fit,power_fit,stddev)

        if   BC == "no_slip":
            plt.plot( x, y,"s",color='darkturquoise',markersize=f_marker*markersize,markerfacecolor='darkturquoise',markeredgecolor='black')
            plt.plot(xx,yy,"-",color='darkturquoise',label=labf,linewidth=f_line*linewidth)
        elif BC == "free_slip":
            if paper == "PEPI-method":
                plt.plot( x, y,"s",color='darkturquoise',markersize=f_marker*markersize,markerfacecolor='darkturquoise',markeredgecolor='black')
                plt.plot(xx,yy,"-",color='darkturquoise',label=labf,linewidth=f_line*linewidth)
            else:
                plt.plot( x, y,"s",color='orange',markersize=f_marker*markersize,markerfacecolor='orange',markeredgecolor='black')
                plt.plot(xx,yy,"-",color='orange',label=labf,linewidth=f_line*linewidth)

    if plot_GL_Iu:

        log_Re_GLIu_min      = np.min(np.log10(np.min(Re_GLIu)))
        log_Re_min_scaling   = min(log_yy_min,log_Re_GLIu_min)

        log_Re_GLIu_max        = np.max(np.log10(np.max(Re_GLIu)))
        log_Re_max_scaling   = max(log_yy_max,log_Re_GLIu_max)

        log_Re_range_scaling = log_Re_max_scaling - log_Re_min_scaling

        log_Re_max_scaling  += 0.01*log_Re_range_scaling
        log_Re_min_scaling  += 0.25*log_Re_range_scaling

        if what == "Re":
            plt.plot(Ra_GL,Re_GLIu,"-",color="gold",linewidth=linewidth,label=r'$Nu_{GL^{I_u}} \ = 0.0.039 Ra^{1/2} Pr^{-5/6}$')

        log_Nu_GL_min        = np.min(np.log10(np.min(Nu_GL)))
        log_Nu_min_scaling   = min(log_yy_min,log_Nu_GL_min)

        log_Nu_GL_max        = np.max(np.log10(np.max(Nu_GL)))
        log_Nu_max_scaling   = max(log_yy_max,log_Nu_GL_max)

        log_Nu_range_scaling = log_Nu_max_scaling - log_Nu_min_scaling
        log_Nu_max_scaling  += 0.48*log_Nu_range_scaling
        log_Nu_min_scaling  -= 0.05*log_Nu_range_scaling

        if what == "Nu":
            plt.plot(Ra_GL,Nu_GL3u,"-",color="orange",linewidth=linewidth,label=r'$Nu_{GL^{III_u}} = Ra^{1/4} Pr^{1/8}\left(0.24+3.3\times 10^{-3}Ra^{5/28}Pr^{-15/56}\right)$')
            plt.plot(Ra_GL,Nu_GL,"-",color="lightgrey",linewidth=linewidth,label=r'$Nu_{GL} \ \ \ = 0.27\ \ \, Ra^{1/4} + 0.038\ \ Ra^{1/3}$')
            plt.plot(Ra_GL,Nu_GLIu,"-",color="gold",linewidth=linewidth,label=r'$Nu_{GL^{I_u}} = 0.33 Ra^{1/4} Pr^{-1/12}$')

    plt.rcParams['xtick.labelsize'] = fontsize_axis_tics
    plt.rcParams['ytick.labelsize'] = fontsize_axis_tics
    plt.tick_params(axis='x',labelsize=fontsize_axis_labels)
    plt.tick_params(axis='y',labelsize=fontsize_axis_labels)

    plt.xlabel(xlab, fontsize=fontsize_axis_labels)
    if what == "Vrms":
        plt.ylabel(r'$V_{rms}$', fontsize=fontsize_axis_labels)
    else:
        plt.ylabel(ylab, fontsize=fontsize_axis_labels)
    plt.xscale("log")
    plt.yscale("log")

#================================================================================
# Main
#================================================================================

#--------------------------------------------------------------------------------
# Loop over all suites of runs to plot
#--------------------------------------------------------------------------------

for iwhat_suite in range(len(whats_suite)):

 what_suite = whats_suite[iwhat_suite]

 for iwhat_version in range(len(whats_version)):

  suite_version = whats_version[iwhat_version]

  if suite_version == "v6": Ras = ["Ra4"]*6

  Nu_bar    = np.zeros((len(BCs),nruns_max))
  Re_bar    = np.zeros((len(BCs),nruns_max))
  Vrms_bar  = np.zeros((len(BCs),nruns_max))
  Ma_bar    = np.zeros((len(BCs),nruns_max))
  Tbar_bar  = np.zeros((len(BCs),nruns_max))
  TBLW_bar  = np.zeros((len(BCs),nruns_max))

  nt_int       = np.zeros((len(BCs),nruns_max),dtype=int)
  nt_run       = np.zeros((len(BCs),nruns_max),dtype=int)
  nt_run_int   = np.zeros((len(BCs),nruns_max),dtype=int)

  for iBC in range(len(BCs)):

    BC = BCs[iBC]

    print("==========================================================================================")
    print("what_suite = %s, version = %s, BCs = %s" %(what_suite,suite_version,BC))
    print("==========================================================================================")

#   Specification of:
#
#   (1) The name to identify the suite, the project containing the runs in the suite, the suite title, and the output filename base
#   (2) Run names in the suite, their colors on time series plots, and whether to label them on time series plots
#   (3) Start and end of normalized time of averaging window for mean Nu, Vrms etc, and for scaling analysis

    if what_suite == "Pr=1000":

        suite_title     = what_suite + ", " + suite_version + resolution + ", " + BC
        suite_extension = suite_version

        if paper != "":
            if suite_extension == "":
                output_filename_base = paper + "_" + what_suite + resolution + "_" + BC
            else:
                output_filename_base = paper + "_" + what_suite + "_"  + suite_version + resolution + "_" + BC
        else          :
            output_filename_base = what_suite + "_"  + suite_version + resolution + "_" + BC

        if test: output_filename_base = output_filename_base + resolution

        runs          = [""]*len(Ras_suite)

        for i in range(len(Ras_suite)): runs[i] = "Ra" + str(Ras_suite[i]) + "_" + BC + "_" + what_suite

        if test:
            for i in range(len(Ras_suite)):
                runs[i] = runs[i] + resolution

        output_filename_base = output_filename_base + "_" + str(averaging_window_start) + "-" + str(averaging_window_end)

#       Runs at half power of 10 increments

        if   len(runs) == 1:  cols = ["blue"]
        elif len(runs) == 2:  cols = ["blue","red"]
        elif len(runs) == 3:  cols = ["blue","gold","red"]
        elif len(runs) == 4:  cols = ["blue","aqua","gold","red"]
        elif len(runs) == 5:  cols = ["darkblue","blue","aqua","gold","red"]
        elif len(runs) == 6:  cols = ["darkblue","blue","aqua","gold","orange","red"]
        elif len(runs) == 7:  cols = ["darkblue","blue","aqua","lightgreen","gold","orange","red"]
        elif len(runs) == 8:  cols = ["darkblue","blue","aqua","lightgreen","greenyellow","gold","orange","red"]
        elif len(runs) == 9:  cols = ["darkblue","blue","dodgerblue","aqua","lightgreen","greenyellow","gold","orange","red"]
        elif len(runs) == 10: cols = ["darkblue","blue","dodgerblue","darkturquoise","aqua","lightgreen","greenyellow","gold","orange","red"]
        elif len(runs) == 11: cols = ["darkblue","blue","dodgerblue","darkturquoise","aqua","lightgreen","greenyellow","gold","orange","darkorange","red"]
        elif len(runs) == 12: cols = ["darkblue","blue","dodgerblue","darkturquoise","aqua","lightgreen","greenyellow","gold","orange","darkorange","red", "darkred"]
        elif len(runs) == 13: cols = [               "darkblue","blue","dodgerblue","deepskyblue","darkturquoise","aqua","lightgreen","greenyellow","gold","orange","darkorange","red","darkred"]
        elif len(runs) == 14: cols = ["midnightblue","darkblue","blue","dodgerblue","deepskyblue","darkturquoise","aqua","lightgreen","greenyellow","gold","orange","darkorange","red","darkred"]
        elif len(runs) == 15: cols = ["midnightblue","darkblue","blue","dodgerblue","deepskyblue","darkturquoise","aqua","lightgreen","greenyellow","gold","orange","darkorange","red","darkred","black"]
        elif len(runs) == 16: cols = ["midnightblue","mediumblue","darkblue","blue","dodgerblue","deepskyblue","darkturquoise","aqua","lightgreen","greenyellow","gold","orange","darkorange","red","darkred","black"]
        elif len(runs) == 17: cols = ["midnightblue","mediumblue","darkblue","blue","dodgerblue","deepskyblue","darkturquoise","aqua","lightgreen","mediumspringgreen","greenyellow","gold","orange","darkorange","red","darkred","black"]
        elif len(runs) == 18: cols = ["midnightblue","mediumblue","darkblue","blue","dodgerblue","deepskyblue","darkturquoise","aqua","lightgreen","mediumspringgreen","green","greenyellow","gold","orange","darkorange","red","darkred","black"]
        elif len(runs) == 19: cols = ["midnightblue","mediumblue","darkblue","blue","dodgerblue","deepskyblue","darkturquoise","aqua","lightgreen","mediumspringgreen","green","greenyellow","yellow","gold","orange","darkorange","red","darkred","black"]
        elif len(runs) == 20: cols = ["midnightblue","mediumblue","darkblue","blue","dodgerblue","deepskyblue","darkturquoise","aqua","lightgreen","mediumspringgreen","green","lime","greenyellow","yellow","gold","orange","darkorange","red","darkred","black"]
        elif len(runs) == 21: cols = ["midnightblue","mediumblue","darkblue","blue","dodgerblue","deepskyblue","darkturquoise","aqua","lightgreen","mediumspringgreen","green","lime","limegreen","greenyellow","yellow","gold","orange","darkorange","red","darkred","black"]
        elif len(runs) == 22: cols = ["black","midnightblue","mediumblue","darkblue","blue","dodgerblue","deepskyblue","darkturquoise","aqua","lightgreen","mediumspringgreen","green","lime","limegreen","greenyellow","yellow","gold","orange","darkorange","red","darkred","black"]
        elif len(runs) == 22: cols = ["midnightblue","mediumblue","darkblue","blue","dodgerblue","deepskyblue","darkturquoise","aqua","lightgreen","mediumspringgreen","lime","limegreen","green","lawngreen","greenyellow","yellow","gold","orange","darkorange","red","darkred","black"]
        elif len(runs) == 23: cols = ["black","midnightblue","mediumblue","darkblue","blue","dodgerblue","deepskyblue","darkturquoise","aqua","lightgreen","mediumspringgreen","lime","limegreen","green","lawngreen","greenyellow","yellow","gold","orange","darkorange","red","darkred","black"]
        elif len(runs) == 24: cols = ["black","midnightblue","mediumblue","darkblue","blue","dodgerblue","deepskyblue","darkturquoise","aqua","lightgreen","mediumspringgreen","lime","limegreen","green","lawngreen","greenyellow","yellow","gold","orange","darkorange","red","darkred","black"]
        else:
            print("ERROR: len(runs) > 23 exceeds defined colours for the plots! Need to add more colors to the cols[] array.")
            exit()

        plot_legends = [True]*len(runs)
        if len(runs) > 11:
            for i in range(1,len(runs),2): plot_legends[i] = False

    elif what_suite == "Pr=1_free_slip":

        project = "benchmarks_annulus"
        suite_title   = "Pr=1 (free slip BC)"
        output_filename_base = "free-slip_Pr=1"
        runs          = ["Ra9_free_slip_Pr=1"  , "Ra10_free_slip_Pr=1"  , "Ra11_free_slip_Pr=1", "Ra12_free_slip_Pr=1"]
        cols          = ["blue"                , "aqua"                 , "gold"               , "red"                ]
        plot_legends = [True]*len(runs)

#   Set the input and output directory paths and create the output directory path if necessary

    in_dir_path          = "./projects/" + project + "/"
    out_dir_path         = in_dir_path   + "bmk/"
    if not os.path.exists(out_dir_path):  os.mkdir(out_dir_path)

#--------------------------------------------------------------------------------
#   Loop over variables to plot as time series
#--------------------------------------------------------------------------------

    print("Runs  = ",runs)
    print("==========================================================================================")

    Ras          = np.zeros(len(runs))
    log_Ra       = np.zeros(len(runs))
    Prs          = np.zeros(len(runs))
    log_Pr       = np.zeros(len(runs))
    As           = np.zeros(len(runs))
    log_A        = np.zeros(len(runs))

    timescale_nt     = ["-"]*len(runs)
    timescale_nt_run = ["-----"]*len(runs)

    nt_in         = np.zeros(len(runs) ,dtype=int)
    dt_p_spec_run = np.zeros(len(runs))
    dt_p_run      = np.zeros(len(runs))

    print("Making time series plots  ... ")
    print("==========================================================================================")

    for iwhat in range(len(whats_tseries)):

        what = whats_tseries[iwhat]

#       Initialize log_10 of dimensionless parameters to calculate the min and max of these for all runs

        log_Ra_min_plot   = 100; log_Ra_max_plot   = -100
        log_Pr_min_plot   = 100; log_Pr_max_plot   = -100
        log_A_min_plot    = 100; log_A_max_plot    = -100

        log_Nu_min_plot   = 100; log_Nu_max_plot   = -100
        log_Vrms_min_plot = 100; log_Vrms_max_plot = -100
        log_Re_min_plot   = 100; log_Re_max_plot   = -100
        log_Ma_min_plot   = 100; log_Ma_max_plot   = -100
        Tbar_min_plot     = 100; Tbar_max_plot     = -100
        TBLW_min_plot     = 100; TBLW_max_plot     = -100

#       Loop over all runs to get information on the suite of runs

        for i in range(len(runs)-1,-1,-1):

#           Set the input file for the run and read it to get the run parameters

            run             = runs[i]
            in_run_file     = in_dir_path + run + ".in"
            get_parameters(in_run_file)
            from get_parameters import *

#           Get R_p and L_p from descriptive names if needed

            R_p, L_p = get_R_and_L(R_p,L_p,False,0)
            if geometry == "cartesian":
                L_l  = nz-1
                dx_p = L_p/L_l
                As[i] = (nx-1)/(nz-1)
            elif geometry == "annulus":
                radii, thetas, nz, L_l, dx, As[i] = get_annular_lattice(R_p,L_p,nx,theta1,theta2,BC_sides,False,0)
                dx_p = L_p/L_l
            error = False

#           Get dt_p for the run and initialize some values needed by get_lattice_viscosity()

            if state == "solid": pturb = "none"
            else:                pturb = "default"
            nu_l_min = 0.0001
            Ma_sim   = "default"
            u_p_geo  = "default"
            gamma_T  = 0
            dt_p_spec_run = dt_p
            optimize_nu = True
            Ra_p, Pr_p, eta_p, nu_p, ka_p, Ra_l, Pr_l, g_l, nu_l, ka_l, dx_p, dt_p_run[i], pturb = get_lattice_viscosity(state,pturb,Ra,Pr,g_p,beta_p,dT_p,L_p,eta_p,rho_p,ka_p,g_l,L_l,nu_l,optimize_nu,nu_l_min,dt_p,Ma_sim,u_p_geo,gamma_T,rheology,geometry,False,0)

#           Get timesteps, nt_int and the nt timescale from the input nt specification

            nt_in[i], nt_int[iBC,i]   , timescale_nt[i]    = get_time_value_text_int_timescale_from_time_spec("nt"      ,nt      ,dt_p_run[i])


#           Read each time series, get nt, and calculate the min and max values

            nt_max   = 10000
            Nu       = np.zeros(nt_max)
            Re       = np.zeros(nt_max)
            Vrms     = np.zeros(nt_max)
            Ma       = np.zeros(nt_max)
            Tbar     = np.zeros(nt_max)
            TBLW     = np.zeros(nt_max)

            it = 0
            in_bmk_file = in_dir_path + run + ".bmk"
            with open(in_bmk_file,"r") as file:

                for line in file:
                    cleaned_line = line.strip()
                    words        = cleaned_line.split()

                    if it == 0:

                        get_indices(words)

                        log_Ra[i] = float(words     [i_Ra])
                        Ras   [i] = float(10**log_Ra[i   ])

                        log_Pr[i] = float(words     [i_Pr])
                        Prs   [i] = float(10**log_Pr[i   ])

                        As    [i] = float(words     [i_A ])
                        log_A [i] = np.log10(As[i])

                    if it > nt_max:
                        print("    ERROR: it = %d > %d = nt_in, len(Nu) = %d" %(it,nt_in[i],len(Nu)))
                        print("           file       : %s" %(in_bmk_file))
                        print("           line number: %d" %(it+1))
                        print("           line       : %s" %(words))
                        exit()

                    nt_run_int  [iBC,i ] = int  (words[i_t_p_int      ])
                    timescale_nt_run[i ] = str  (words[i_timescale_t_p])
                    Nu              [it] = float(words [i_Nu          ])
                    Re              [it] = float(words[i_Re           ])
                    Vrms            [it] = float(words[i_Vrms         ])
                    Tbar            [it] = float(words[i_Tbar         ])
                    Ma              [it] = float(words[i_Ma           ])
                    TBLW            [it] = np.abs(float(words[i_TBLW  ]))

                    it           += 1

            nt_run[iBC,i] = it

            Nu    = Nu  [:nt_run[iBC,i]]
            Re    = Re  [:nt_run[iBC,i]]
            Vrms  = Vrms[:nt_run[iBC,i]]
            Ma    = Ma  [:nt_run[iBC,i]]
            Tbar  = Tbar[:nt_run[iBC,i]]
            TBLW  = TBLW[:nt_run[iBC,i]]

            it0 = int((nt_run[iBC,i]-1)*averaging_window_start+0.5)  # Start of averaging window
            it0 = max(0        ,it0)
            it1 = int((nt_run[iBC,i]-1)*averaging_window_end  +0.5)  # Start of averaging window
            it1 = min(nt_run[iBC,i],it1)

#           Update min and max log range of dimensionless parameters of runs, Ra, Pr and A

            log_Ra_min_plot = min(np.min(log_Ra),log_Ra_min_plot)
            log_Ra_max_plot = max(np.max(log_Ra),log_Ra_max_plot)
    
            log_Pr_min_plot = min(np.min(log_Pr),log_Pr_min_plot)
            log_Pr_max_plot = max(np.max(log_Pr),log_Pr_max_plot)

            log_A_min_plot  = min(np.min(log_A ),log_A_min_plot )
            log_A_max_plot  = max(np.max(log_A ),log_A_max_plot )

#           Update min and max log ranges of Nu, Re and Vrms for the plot to min and max of all series

            log_Nu            = np.log10(Nu[nt_run[iBC,i]//5:])

            log_Nu_min        = np.min(log_Nu  )
            log_Nu_min_plot   = min(log_Nu_min  ,log_Nu_min_plot  )

            log_Nu_max        = np.max(log_Nu)
            log_Nu_max_plot   = max(log_Nu_max  ,log_Nu_max_plot  )

            log_Re            = np.log10(Re)

            log_Re_min        = np.min(log_Re  )
            log_Re_min_plot   = min(log_Re_min  ,log_Re_min_plot  )

            log_Re_max        = np.max(log_Re  )
            log_Re_max_plot   = max(log_Re_max,log_Re_max_plot)

            log_Vrms          = np.log10(Vrms)

            log_Vrms_min      = np.min(log_Vrms[nt_run[iBC,i]//20:])
            log_Vrms_min_plot = min(log_Vrms_min,log_Vrms_min_plot)

            log_Vrms_max      = np.max(log_Vrms[nt_run[iBC,i]//20:])
            log_Vrms_max_plot = max(log_Vrms_max,log_Vrms_max_plot)

            log_Ma            = np.log10(Ma)
            log_Ma_min        = np.min(log_Ma  [nt_run[iBC,i]//2:])
            log_Ma_min_plot   = min(log_Ma_min,log_Ma_min_plot)

            log_Ma_max        = np.max(log_Ma  [nt_run[iBC,i]//2:])
            log_Ma_max_plot   = max(log_Ma_max,log_Ma_max_plot)

            Tbar_min          = np.min(Tbar    [nt_run[iBC,i]//50:])
            Tbar_min_plot     = min(Tbar_min,Tbar_min_plot)

            Tbar_max          = np.max(Tbar    [nt_run[iBC,i]//50:])
            Tbar_max_plot     = max(Tbar_max,Tbar_max_plot)

            TBLW_min          = np.min(TBLW    [nt_run[iBC,i]//2:])
            TBLW_min_plot     = min(TBLW_min,TBLW_min_plot)

            TBLW_max          = np.max(TBLW    [nt_run[iBC,i]//2:])
            TBLW_max_plot     = max(TBLW_max,TBLW_max_plot)

#       Set the figure size and time series plot title

        fig = plt.figure(figsize=(xsize,ysize))  ; plt.subplots_adjust(left=left,bottom=bottom,top=top,right=right)

        if what == "Tbar": title = r'$\overline{T}(t)$'         + ", %s" %(suite_title)
        else             : title = r'$%s          (t)$' %(what) + ", %s" %(suite_title)

        if paper == "GJI-benchmark" or paper == "PEPI-method": # PEPI method paper and GJI benchmark paper titles

            if   what == "Nu"  : title = "Nusselt number versus time"
            elif what == "Vrms": title = r'$V_{rms}$' + " versus time"
            elif what == "Re"  : title = "Reynolds number versus time"
            elif what == "Ma"  : title = "Mach number versus time"
            elif what == "TBLW": title = "TBL width versus time"
            elif what == "Tbar": title = r'$\bar T$' + " versus time"
            if   BC == "free_slip":
                title = title + " (free-slip)"
            elif BC == "no_slip":
                title = title + " (zero-slip)"
            if suite_version != "": title = title + " (%s)" %(suite_version)

        plt.title(title,fontsize=fontsize_title)

#       Set axis tick sizes

        plt.xlabel(r'Normalized time',fontsize=fontsize_axis_labels)

        plt.rcParams['xtick.labelsize'] = fontsize_axis_tics
        plt.rcParams['ytick.labelsize'] = fontsize_axis_tics

        plt.tick_params(axis='x',labelsize=fontsize_axis_labels)
        plt.tick_params(axis='y',labelsize=fontsize_axis_labels)

#       Set y limits of plot

        plt.xlim(0,averaging_window_end)
        if   what == "Nu":
            log_yy_min_plot = log_Nu_min_plot   ; log_yy_max_plot = log_Nu_max_plot
        elif what == "Re":
            log_yy_min_plot = log_Re_min_plot   ; log_yy_max_plot = log_Re_max_plot
        elif what == "Vrms":
            log_yy_min_plot = log_Vrms_min_plot ; log_yy_max_plot = log_Vrms_max_plot
        elif what == "Ma":
            log_yy_min_plot = log_Ma_min_plot   ; log_yy_max_plot = log_Ma_max_plot
        elif what == "Tbar":
            yy_min_plot     = Tbar_min_plot     ; yy_max_plot     = Tbar_max_plot
        elif what == "TBLW":
            log_yy_min_plot = np.log10(TBLW_min_plot) ; log_yy_max_plot = np.log10(TBLW_max_plot)

        if what != "Tbar":
            log_yy_range_plot = log_yy_max_plot-log_yy_min_plot
            if   what == "Vrms": factor = 1.3
            elif what == "Ma"  : factor = 1.22
            elif what == "TBLW": factor = 1.28
            elif what == "Nu"  : factor = 1.15
            else               : factor = 1
            if len(runs) <= 6: 
                log_yy_max_plot = log_yy_min_plot + (1+len(runs)/20)*(factor*log_yy_range_plot)
            else:
                log_yy_max_plot = log_yy_min_plot + factor*log_yy_range_plot
            log_yy_range_plot = log_yy_max_plot-log_yy_min_plot
            log_yy_min_plot -= 0.05*log_yy_range_plot
            yy_min_plot  = 10**log_yy_min_plot  ; yy_max_plot     = 10**log_yy_max_plot
        elif what == "Tbar":
            yy_range_plot = yy_max_plot-yy_min_plot
            factor = 1.2
            yy_max_plot = yy_min_plot + (1+len(runs)/20)*(factor*yy_range_plot)
            yy_range_plot = yy_max_plot-yy_min_plot
            yy_min_plot -= 0.05*yy_range_plot

#--------------------------------------------------------------------------------
#       Loop again over runs, select each time series, calculate the averages over
#       the analysis window for the scaling plot, and plot each time series
#--------------------------------------------------------------------------------

        for i in range(len(runs)-1,-1,-1):

            run             = runs[i]

            in_bmk_file     = in_dir_path + run + ".bmk"
            in_run_file     = in_dir_path + run + ".in"

            it0 = int((nt_run[iBC,i]-1)*averaging_window_start+0.5)  # Start of averaging window
            it0 = max(0        ,it0)
            it1 = int((nt_run[iBC,i]-1)*averaging_window_end  +0.5)  # Start of averaging window
            it1 = min(nt_run[iBC,i],it1)

            Nu   = np.zeros(nt_run[iBC,i])
            Re   = np.zeros(nt_run[iBC,i])
            Vrms = np.zeros(nt_run[iBC,i])
            Ma   = np.zeros(nt_run[iBC,i])
            Tbar = np.zeros(nt_run[iBC,i])
            TBLW = np.zeros(nt_run[iBC,i])

            it    = 0
            n_sum = 0
            with open(in_bmk_file,"r") as file:

                for line in file:
                    cleaned_line = line.strip()
                    words        = cleaned_line.split()

                    if it == 0: get_indices(words)

#                   Save the time series of Nu, Re and Vrms and calculate log_10 of
#                   each series along with their averages over the analysis window

                    Nu      [it]  = float(words[i_Nu       ])
                    Re      [it]  = float(words[i_Re       ])
                    Vrms    [it]  = float(words[i_Vrms     ])
                    Ma      [it]  = float(words[i_Ma       ])
                    Tbar    [it]  = float(words[i_Tbar     ])
                    TBLW    [it]  = np.abs(float(words[i_TBLW     ]))

                    if it >= it0 and it <= it1:
                        Nu_bar   [iBC,i] += Nu  [it]
                        Re_bar   [iBC,i] += Re  [it]
                        Vrms_bar [iBC,i] += Vrms[it]
                        Ma_bar   [iBC,i] += Ma  [it]
                        Tbar_bar [iBC,i] += Tbar[it]
                        TBLW_bar [iBC,i] += TBLW[it]
                        n_sum       += 1

                    it += 1

#               Calculate the mean of Nu, Re and Vrms over the analysis window

            if n_sum > 0:
                Nu_bar    [iBC,i] /= n_sum
                Re_bar    [iBC,i] /= n_sum
                Vrms_bar  [iBC,i] /= n_sum
                Ma_bar    [iBC,i] /= n_sum
                Tbar_bar  [iBC,i] /= n_sum
                TBLW_bar  [iBC,i] /= n_sum

            Nu    = Nu   [:nt_run[iBC,i]]
            Re    = Re   [:nt_run[iBC,i]]
            Vrms  = Vrms [:nt_run[iBC,i]]
            Ma    = Ma   [:nt_run[iBC,i]]
            Tbar  = Tbar [:nt_run[iBC,i]]
            TBLW  = TBLW [:nt_run[iBC,i]]

#           Set log y scale for Nu, Re and Vrms and what time series will be plotted, and plot the y label

            if what != "Tbar": plt.yscale("log")

            if n_sum == 0:
                lab_bar = r'$\overline{Nu}$' + "unknown"
            else:
                if   what == "Nu"  :
                    f = Nu
                    plt.ylabel(r'$Nu$'     ,fontsize=fontsize_axis_labels)
#                    if   Nu_bar[iBC,i] <   10: lab_bar = r'$\overline{Nu} = %5.3f$' %(Nu_bar[iBC,i])
#                    elif Nu_bar[iBC,i] <  100: lab_bar = r'$\overline{Nu} = %6.2f$' %(Nu_bar[iBC,i])
#                    elif Nu_bar[iBC,i] < 1000: lab_bar = r'$\overline{Nu} = %7.1f$' %(Nu_bar[iBC,i])
#                    else                 : lab_bar = r'$\overline{Nu} = %9.4g$' %(Nu_bar[iBC,i])
                    lab_bar = r'$\overline{Nu} = %.4g$' %(Nu_bar[iBC,i])
                elif what == "Re"  :
                    f = Re
                    plt.ylabel(r'$Re$'     ,fontsize=fontsize_axis_labels)
                    lab_bar = r'$\overline{Re} = %.4g$' %(Re_bar[iBC,i])
                elif what == "Vrms":
                    f = Vrms
                    plt.ylabel(r'$V_{rms}$',fontsize=fontsize_axis_labels)
                    lab_bar = r'$\overline{V}_{rms} = %.4g$' %(Vrms_bar[iBC,i])
                elif what == "Ma"  :
                    f = Ma
                    plt.ylabel(r'$Ma$',fontsize=fontsize_axis_labels)
                    lab_bar = r'$\overline{Ma} = %6.4g$' %(Ma_bar[iBC,i])
                elif what == "Tbar":
                    f = Tbar
                    plt.ylabel(r'$\bar T$'     ,fontsize=fontsize_axis_labels)
                    lab_bar = r'$\overline{\bar T} = %5.3f$' %(Tbar_bar[iBC,i])
                elif what == "TBLW":
                    f = TBLW
                    plt.ylabel(r'$TBLW$'     ,fontsize=fontsize_axis_labels)
                    lab_bar = r'$\overline{TBLW} = %5.4g$' %(TBLW_bar[iBC,i])

#           Physical time of run as a string

            t_p_run = timevalue_to_secs(nt_run_int[iBC,i],timescale_nt_run[i],0)
            lab_t_p_run = time_to_physical_time_text(t_p_run,dt_p_run[i])

            if t_p_run/sec_gy   >= 1: lab_t_p_run = "%4.1f Gy" %(t_p_run/sec_gy)
            elif t_p_run/sec_my >= 1: lab_t_p_run = "%4.1f My" %(t_p_run/sec_my)

#           Set the y limits of the plot

            if what != "Tbar" and what != "TBLW": plt.ylim(10**(log_yy_min_plot),10**(log_yy_max_plot))
            else                                : plt.ylim(         yy_min_plot ,         yy_max_plot )

#           Calculate the factor to convert physical to normalized time and plot the time series

            t_normalized = np.linspace(0,1,nt_run[iBC,i])
            if plot_legends[i]: plt.plot(t_normalized,f,"-",color=cols[i],linewidth=linewidth,label=r'$Ra=10^{%4.2f}$' %(log_Ra[i])+", " + lab_bar + " (" + "%s" %(lab_t_p_run) + ")")
            else              : plt.plot(t_normalized,f,"-",color=cols[i],linewidth=linewidth)

#           Plot the start of the averaging window and legend

            if i == 0:
                plt.axvline(averaging_window_start,linestyle="-",color="black",linewidth=1)
                if averaging_window_end < 1:
                    plt.axvline(averaging_window_end  ,linestyle="--",color="black",linewidth=1)

            plt.legend(fontsize=fontsize_legend_tseries,loc="upper left")

#       Save the time series plot

        if iwhat == 0:
            for i in range(len(Ras)):
                print("Ra = 10^%3.1f, Nu_bar = %8.4g, Vrms_bar = %8.4g, Tbar = %5.3f, Ma_bar = %8.3g, TBLW = %8.3g" %(np.log10(Ras[i]),Nu_bar[iBC,i],Vrms_bar[iBC,i],Tbar_bar[iBC,i],Ma_bar[iBC,i],TBLW_bar[iBC,i]))
            print("------------------------------------------------------------------------------------------")
        filename      = out_dir_path + output_filename_base + "_tseries_%s" %(what)
        full_filename = filename + '.png'
        try: os.remove(full_filename)
        except OSError: pass
        plt.savefig(full_filename)
        plt.close()
        print(full_filename)

#--------------------------------------------------------------------------------
#   Loop for analysis of scaling relations and to generate scaling plots
#--------------------------------------------------------------------------------

#    Calculate Grossman-Lohse scaling for Nu to Re

    n_GL = 101
    dlogRa_GL = (log_Ra_max_plot-log_Ra_min_plot)/(n_GL-1)
    Nu_GL  = np.zeros(n_GL)
    Nu_GL3u = np.zeros(n_GL)
    Ra_GL  = np.zeros(n_GL)
    Nu_GLIu  = np.zeros(n_GL)
    Re_GLIu  = np.zeros(n_GL)
    for i_GL in range(n_GL):
        Ra_GL[i_GL] = 10**(log_Ra_min_plot+i_GL*dlogRa_GL)
        Nu_GL[i_GL] = 0.27*Ra_GL[i_GL]**(1/4)+0.038*Ra_GL[i_GL]**(1/3) # Calculate Nu_GL
        Nu_GL3u[i_GL] = Ra_GL[i_GL]**(1/4)*Prs[0]**(1/8)* ( 0.24 + 3.3e-3*Ra_GL[i_GL]**(5/28)*Prs[0]**(-15/56) )
        Nu_GLIu[i_GL] = 0.33*Ra_GL[i_GL]**(1/4)*Prs[0]**(-1/12)
        Re_GLIu[i_GL] = 0.039*Ra_GL[i_GL]**(1/2)*Prs[0]**(-5/6)

if len(whats_scaling) > 0:
    print("==========================================================================================")
    print("Making scaling plots  ... ")
    print("==========================================================================================")

for iwhat_suite in range(len(whats_suite)):

 what_suite = whats_suite[iwhat_suite]

 for iwhat_version in range(len(whats_version)):

   suite_version = whats_version[iwhat_version]

   if suite_version == "v6": Ras = ["Ra4"]*6

   if paper != "":
       if paper == "GJI-benchmark":
           output_filename_base = paper + "_" + what_suite
           if suite_version != "": output_filename_base + "_"  + suite_version + resolution
       elif suite_version == "":
           output_filename_base = paper + "_" + what_suite + resolution + "_" + BC
       else:
           output_filename_base = paper + "_" + what_suite + "_"  + suite_version + resolution + "_" + BC
   else          :
       output_filename_base = what_suite + "_"  + suite_version + resolution + "_" + BC
  
#   for j in range(len(BCs)): output_filename_base = output_filename_base + "_" + BCs[j]

   for iwhat_scaling in range(len(whats_scaling)):

     what = whats_scaling[iwhat_scaling]
     filename  = output_filename_base  + "_scaling_%s" %(what)
     filename  = out_dir_path + filename + '.png'
     title = "%s, " %(suite_title) + r'$%s - Ra$' %(what)  + " scaling"
     if suite_version != "": title = title + " (%s)" %(suite_version)
     try: os.remove(filename)
     except OSError: pass

     if paper == "GJI-benchmark" or paper == "PEPI-method": # PEPI method paper and GJI benchmark paper titles
         if   what == "Nu"  : title = "Nusselt number scaling"
         elif what == "Vrms": title = "RMS velocity scaling"
         elif what == "Re"  : title = "Reynolds number scaling"
     else:
         title = title + " (%s)" %(suite_version)

     fig = plt.figure(figsize=(xsize,ysize)) ; plt.subplots_adjust(left=left_sc,bottom=bottom_sc,top=top_sc,right=right_sc)
     plt.title(title,fontsize=fontsize_title)

     log_xx_min  = 100 ; log_xx_max  = -100
     log_yy_min  = 100 ; log_yy_max  = -100

     for iBC in range(len(BCs)):

       BC = BCs[iBC]

       runs          = [""]*len(Ras_suite)

       for i in range(len(Ras_suite)): runs[i] = "Ra" + str(Ras_suite[i]) + "_" + BC + "_" + what_suite

       if test:
           for i in range(len(Ras_suite)):
               runs[i] = runs[i] + resolution

       print("-------------- ",BC," ---------------")
       for icode in range(len(codes)):

        code = codes[icode]

        if code == "GTXLBM":
          f_bar  = np.zeros(len(runs))
          x_data = np.zeros(len(runs))
          y_data = np.zeros(len(runs))
          y_fit  = np.zeros(len(runs))
          n_data = 0
          for i in range(len(runs)):
            if nt_run[iBC,i] > 1:
                x_data[n_data] = log_Ra[i]
                if   what == "Re" and Re_bar[iBC,i] > 0:
                    y_data[n_data] = np.log10(Re_bar  [iBC,i])
                elif what == "Nu" and Nu_bar[iBC,i] > 0:
                    y_data[n_data] = np.log10(Nu_bar  [iBC,i])
                elif what == "Vrms" and Vrms_bar[iBC,i] > 0:
                    y_data[n_data] = np.log10(Vrms_bar[iBC,i])
                elif what == "Pr"   and Prs           [    i] > 0:
                    y_data[n_data] = np.log10(Prs     [    i])
                elif what == "A"    and As            [    i] > 0:
                    y_data[n_data] = np.log10(As      [    i])
                n_data += 1

          x_data = x_data[:n_data]
          y_data = y_data[:n_data]

          c1,c0 = np.polyfit(x_data,y_data, 1)

          y_fit = c0 + c1*x_data

          variance_sum = 0
          for j in range(len(y_data)):
              dvariance = (np.abs(y_data[j]-y_fit[j])/y_fit[j])**2
              variance_sum += dvariance
          variance = variance_sum/len(y_data)
          stddev = np.sqrt(variance)

          const_fit = 10**c0
          power_fit = c1

          plot_scaling(title,filename,10**x_data,10**y_data,10**x_data,10**y_fit,r'$Ra$',r'$%s$' %(what))

        elif code == "ASPECT" and (what == "Nu" or what == "Vrms"):

            if   BC == "no_slip":
                if   what == "Nu"  : y_data = Nus_ASPECT_no_slip
                elif what == "Vrms": y_data = Vrms_ASPECT_no_slip
            elif BC == "free_slip":
                if   what == "Nu"  : y_data = Nus_ASPECT_free_slip
                elif what == "Vrms": y_data = Vrms_ASPECT_free_slip
            x_data = np.log10(Ras_ASPECT)
            y_data = np.log10(y_data)

            c1,c0 = np.polyfit(x_data,y_data, 1)

            y_fit = c0 + c1*x_data

            variance_sum = 0
            for j in range(len(y_data)):
                dvariance = (np.abs(y_data[j]-y_fit[j])/y_fit[j])**2
                variance_sum += dvariance
            variance = variance_sum/len(y_data)
            stddev = np.sqrt(variance)

            const_fit = 10**c0
            power_fit = c1

            plot_scaling(title,filename,10**x_data,10**y_data,10**x_data,10**y_fit,r'$Ra$',r'$%s$' %(what))

        error     = (y_data - y_fit)/y_fit
        rms_error = np.sqrt(np.sum(error**2)/len(y_fit))
        print("%s: %s(Ra) = %8.3g + Ra^(%8.3g +/- %8.3g)" %(code,what,const_fit,power_fit,100*rms_error))
        print("      err(log(Ra)) = data - fit = ",end="")
        for i in range(len(x_data)):
            print("(%3.1f,%4.2f" %(x_data[i],100*error[i]),end="")
            print(" %), ",end="")
        print("gamma = %5.3f +/- %5.3f" %(power_fit,stddev))

       print("==========================================================================================")
     log_xx_range = log_xx_max - log_xx_min
     log_yy_range = log_yy_max - log_yy_min

     log_xx_min -= 0.02*log_xx_range
     log_xx_max += 0.02*log_xx_range

     if paper == "PEPI-method":
         log_yy_max += 0.05*log_yy_range
     else:
         log_yy_max += 0.4*log_yy_range
     if suite_version == "v11":
         log_yy_min -= 0.2*log_yy_range
     else:
         log_yy_min -= 0.05*log_yy_range

     plt.xlim(10**log_xx_min,10**log_xx_max)
     plt.ylim(10**log_yy_min,10**log_yy_max)

     plt.legend(fontsize=fontsize_legend_scaling,loc="upper left")

     plt.savefig(filename)
     plt.close()
     print(filename)
     print("==========================================================================================")

