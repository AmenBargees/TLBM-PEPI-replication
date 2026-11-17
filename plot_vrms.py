#================================================================================
# Set suite versions, suites, BCs, and what time series
# and scaling plots to make, and set figsize below here
#================================================================================

import numpy as np

TBLW_data = True

whats_tseries = ["Vrms"]
whats_suite   = ["van-Keken"]
BCs           = [""]
project       = "RT"
display_legend = True
runs = ["van-keken-2x-optimize-nu"]
whats_version = ["-2x-dt=10yrs"]
#runs = ["van-keken-test-optimize-nu","van-keken-test-optimize-nu-0.1"]
#whats_version = ["-test"]
#runs = ["van-keken-2x"]
#whats_version = ["-no-optimize_nu-2x"]
tscale = np.zeros(len(runs))
nruns_max     = len(runs)

# Plot fontsizes

figsize                    = "small"

if figsize == "small":
    left    = 0.13 ; bottom    = 0.12 ; top    = 0.93 ; right    = 0.975
    left_sc = 0.13 ; bottom_sc = 0.12 ; top_sc = 0.93 ; right_sc = 0.985
    linewidth               = 2.5
    markersize              = 6
    xsize                   = 7
    ysize                   = 5
    fontsize_title          = 14
    fontsize_axis_labels    = 14
    fontsize_axis_tics      = 10
    fontsize_legend_tseries = 10
    fontsize_legend_scaling = 12
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
i_Nu           = 7
i_Vrms     = 8
i_Re       = 9
i_TBLW     = 10
i_Nu_u     = 11
i_Nu_d     = 12

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
# from utilities import get_annular_lattice
from utilities import get_R_and_L
from utilities import get_time_value_text_int_timescale_from_time_spec
from utilities import timesteps_to_time_text
from utilities import timetext
from utilities import timevalue
from utilities import timevalue_to_secs
from utilities import timesteps_to_time_text
from utilities import time_to_physical_time_text
from utilities import timevalue_to_timesteps

import get_parameters 
from get_parameters import get_parameters

#================================================================================
# Main
#================================================================================

#--------------------------------------------------------------------------------
# Loop over all suites of runs to plot
#--------------------------------------------------------------------------------

for iwhat_suite in range(len(whats_suite)):
 what_suite = whats_suite[iwhat_suite]

 for iwhat_version in range(len(whats_version)):
  what_version = whats_version[iwhat_version]

  output_filename_base = what_suite + what_version

  nt_int       = np.zeros((len(BCs),nruns_max),dtype=int)
  nt_run       = np.zeros((len(BCs),nruns_max),dtype=int)
  nt_run_int   = np.zeros((len(BCs),nruns_max),dtype=int)

  for iBC in range(len(BCs)):

    print("==========================================================================================")
    print("what_suite = %s" %(what_suite))
    print("==========================================================================================")

    plot_legends = [display_legend]*len(runs)
    cols = ["red","blue"]

#   Set the input and output directory paths and create the output directory path if necessary

    in_dir_path          = "./projects/" + project + "/"
    out_dir_path         = in_dir_path   + "bmk/"
    if not os.path.exists(out_dir_path):  os.mkdir(out_dir_path)

#--------------------------------------------------------------------------------
#   Loop over variables to plot as time series
#--------------------------------------------------------------------------------

    timescale_nt     = ["-"]*len(runs)
    timescale_nt_run = ["-----"]*len(runs)

    nt_in         = np.zeros(len(runs) ,dtype=int)
    dt_p_spec_run = np.zeros(len(runs))
    dt_p_run      = np.zeros(len(runs))

    print("Making time series plots  ... ")
    print("==========================================================================================")

    for iwhat in range(len(whats_tseries)):

        what = whats_tseries[iwhat]

#       Loop over all runs to get information on the suite of runs

        for i in range(len(runs)):

#           Set the input file for the run and read it to get the run parameters

            run             = runs[i]
            in_run_file     = in_dir_path + run + ".in"
            get_parameters(in_run_file)
            from get_parameters import *
            dt_p_spec_run = time_per_timestep
            dt_p = time_per_timestep

#           Get R_p and L_p from descriptive names if needed

            R_p, L_p = get_R_and_L(R_p,L_p)
            if geometry == "cartesian":
                A = (nx-1)/(nz-1)
                L_l  = nz-1
                dx_p = L_p/L_l
            elif geometry == "annulus":
                radii, thetas, nz, L_l, dx, A = get_annular_lattice(R_p,L_p,nx,theta1,theta2,BC_sides,False,0)
                dx_p = L_p/L_l
            error = False

#           Get dt_p for the run and initialize some values needed by get_lattice_viscosity()

            if state == "solid": pturb = "none"
            else:                pturb = "default"
            nu_l_min = 0.0001
            Ma_sim   = "default"
            u_p_geo  = "default"
            gamma_T  = 0

            D = 2 ; rank = 1
            Ra_p, Pr_p, eta_p, nu_p, ka_p, Ra_l, Pr_l, g_l, nu_l, ka_l, dx_p, dt_p_run[i], dt_dt0 = get_lattice_viscosity(state, power_correction,Ra,Pr,g_p,beta_p,dT_p,L_p,eta_p,rho_p,ka_p,Tau_y0_p,D,g_l,L_l,nu_l,nu_l_min,time_per_timestep,gamma_T,rheology,rank=rank)
#           Ra_p, Pr_p, eta_p, nu_p, ka_p, Ra_l, Pr_l, g_l, nu_l, ka_l, dx_p, dt_p_run[i], pturb = get_lattice_viscosity(state,pturb,Ra,Pr,g_p,beta_p,dT_p,L_p,eta_p,rho_p,ka_p,g_l,L_l,nu_l,nu_l_min,dt_p,Ma_sim,u_p_geo,gamma_T,rheology,geometry)

#           Get timesteps, nt_int and the nt timescale from the input nt specification

#            error = False
#            nt_int, timescale_nt, timescale_error = timevalue("nt"  ,nt  ,error,rank)
#            nt_in[i], nt_int[iBC,i]   , timescale_nt[i]    = get_time_value_text_int_timescale_from_time_spec("nt"      ,nt      ,dt_p_run[i])


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
            Vrms /= np.max(Vrms)

#       Set the figure size and time series plot title

        fig = plt.figure(figsize=(xsize,ysize))  ; plt.subplots_adjust(left=left,bottom=bottom,top=top,right=right)

        if   what == "Nu"  : title = "Nusselt number versus time"
        elif what == "Vrms": title = "RMS velocity versus time"
        elif what == "Re"  : title = "Reynolds number versus time"
        elif what == "Ma"  : title = "Mach number versus time"
        elif what == "TBLW": title = "TBL width versus time"
        elif what == "Tbar": title = r'$\bar T$' + " versus time"

        plt.title(title,fontsize=fontsize_title)

#       Set axis tick sizes

        plt.xlabel(r'Normalized time',fontsize=fontsize_axis_labels)

        plt.rcParams['xtick.labelsize'] = fontsize_axis_tics
        plt.rcParams['ytick.labelsize'] = fontsize_axis_tics

        plt.tick_params(axis='x',labelsize=fontsize_axis_labels)
        plt.tick_params(axis='y',labelsize=fontsize_axis_labels)

#       Set y limits of plot

#        plt.xlim(0,1)

#--------------------------------------------------------------------------------
#       Loop again over runs, select each time series, calculate the averages over
#       the analysis window for the scaling plot, and plot each time series
#--------------------------------------------------------------------------------

        for i in range(len(runs)):

            run             = runs[i]

            in_bmk_file     = in_dir_path + run + ".bmk"
            in_run_file     = in_dir_path + run + ".in"

            Nu   = np.zeros(nt_run[iBC,i])
            Re   = np.zeros(nt_run[iBC,i])
            Vrms = np.zeros(nt_run[iBC,i])
            Ma   = np.zeros(nt_run[iBC,i])
            Tbar = np.zeros(nt_run[iBC,i])
            TBLW = np.zeros(nt_run[iBC,i])

            it    = 0
            with open(in_bmk_file,"r") as file:

                for line in file:
                    cleaned_line = line.strip()
                    words        = cleaned_line.split()

#                   Save the time series of Nu, Re and Vrms and calculate log_10 of each series

                    Nu      [it]  = float(words[i_Nu       ])
                    Re      [it]  = float(words[i_Re       ])
                    Vrms    [it]  = float(words[i_Vrms     ])
                    Ma      [it]  = float(words[i_Ma       ])
                    Tbar    [it]  = float(words[i_Tbar     ])
                    TBLW    [it]  = np.abs(float(words[i_TBLW     ]))

                    it += 1

#           Set log y scale for Nu, Re and Vrms and what time series will be plotted, and plot the y label

            if   what == "Nu"  :
                f = Nu
                plt.ylabel(r'$Nu$'     ,fontsize=fontsize_axis_labels)
            elif what == "Re"  :
                f = Re
                plt.ylabel(r'$Re$'     ,fontsize=fontsize_axis_labels)
            elif what == "Vrms":
                Vrms /= np.max(Vrms)
                f = Vrms
                plt.ylabel(r'$V_{rms}$',fontsize=fontsize_axis_labels)
            elif what == "Ma"  :
                f = Ma
                plt.ylabel(r'$Ma_{max}$',fontsize=fontsize_axis_labels)
            elif what == "Tbar":
                f = Tbar
                plt.ylabel(r'$\bar T$'     ,fontsize=fontsize_axis_labels)
            elif what == "TBLW":
                f = TBLW
                plt.ylabel(r'$TBLW$'     ,fontsize=fontsize_axis_labels)

#           Physical time of run as a string

            t_p_run = timevalue_to_secs(nt_run_int[iBC,i],timescale_nt_run[i],0)
            lab_t_p_run = time_to_physical_time_text(t_p_run,dt_p_run[i])

            if t_p_run/sec_gy   >= 1: lab_t_p_run = "%4.1f Gy" %(t_p_run/sec_gy)
            elif t_p_run/sec_my >= 1: lab_t_p_run = "%4.1f My" %(t_p_run/sec_my)

#           Calculate the factor to convert physical to normalized time and plot the time series


            it_maxima = np.zeros(10,dtype=int)
            t_maxima  = np.zeros(10,dtype=float)
            nt = nt_run[iBC,i]
            n_maxima = 0
            t_edge = max(nt//50,10)
            for it in range(t_edge,nt-t_edge):
                if (f[it] >= f[it-1] and f[it] >= f[it+1]) and \
                   (f[it] >= f[it-2] and f[it] >= f[it+2]) and \
                   (f[it] >= f[it-3] and f[it] >= f[it+3]) and \
                   (f[it] >= f[it-4] and f[it] >= f[it+4]) and \
                   (f[it] >= f[it-5] and f[it] >= f[it+5]):
                    if n_maxima == 0 or it/nt > 1.5*t_maxima[0]:
                        it_maxima[n_maxima] = it
                        t_maxima [n_maxima] = it/nt
                        n_maxima += 1
                        print("i = %d, t_maxima[%d] = %5.3f" %(i,n_maxima,t_maxima[n_maxima-1]))
                        if n_maxima == 1:
                            if i == 0: t_max1 = t_maxima[n_maxima-1]
                            tscale[i] = t_maxima[n_maxima-1]/t_max1
                            print("    i = ",i,", t_max1 = ",t_max1,", tscale = ",tscale[i])
                    if n_maxima == 2: break

            scale_factor = nt_run[iBC,i]/nt_run[iBC,0]*tscale[i]
            if i == 0:
                 scale_factor1 = 1
            if i == 1:
                 scale_factor1 = 0.9
            t_normalized = np.linspace(0,scale_factor*scale_factor1,nt_run[iBC,i])
            if plot_legends[i]: plt.plot(t_normalized,f,"-",color=cols[i],linewidth=linewidth,label=what_suite+"(%s)" %(run))
            else              : plt.plot(t_normalized,f,"-",color=cols[i],linewidth=linewidth)

            linestyles = ["--",":"]
            for i_maxima in range(n_maxima):
                print("i = ",i,", maxima = ",i_maxima+1,", linestyle = ",linestyles[i_maxima])
                plt.axvline(t_maxima[i_maxima]/tscale[i],linestyle=linestyles[i_maxima],color=cols[i],linewidth=0.7*linewidth,label=r'$t=%5.3f$' %(t_maxima[i_maxima]/tscale[i]))

#           Plot the start of the averaging window and legend

            plt.legend(fontsize=fontsize_legend_tseries,loc="upper right")

#       Save the time series plot

        if iwhat == 0:
            print("------------------------------------------------------------------------------------------")
        filename      = out_dir_path + output_filename_base + "_tseries_%s" %(what)
        full_filename = filename + '.png'
        try: os.remove(full_filename)
        except OSError: pass
        plt.savefig(full_filename)
        plt.close()
        print(full_filename)

print("------------------------------------------------------------------------------------------")
