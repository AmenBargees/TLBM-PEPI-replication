#===============================================================================
# Generate plot files (png) from data files (npy)
#===============================================================================

import os
from os.path import exists
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
from sys import argv
from datetime import datetime

LaTeX = False
if LaTeX:
    matplotlib.rcParams['text.usetex'] = True

#--------------------------------------------------------------------------------

def plot_graph_line(title,filename,x_i,y_i,xx,yy,xlab,ylab):

    fig = plt.figure(1,figsize=(4.5,3.5))
    plt.subplots_adjust(left=0.19, bottom=0.15, top=0.92, right=0.97)
    plt.title(title,fontsize=12)

#   b = blue
#   k = black
#   w = white

    labd = "Data"
    labf = "Fit"
    col = "blue"
    plt.plot(y_i, x_i,"o",color=col,label=labd,markersize=4)
#    plt.plot(y_i[2:n-2], np.log10(Nu_GL[2:n-2]),"b-",linewidth=1,label="GL theory")
    plt.plot(xx,yy,"-",label=labf + r', $\gamma = $' + "%5.3f" %(gamma),color=col,linewidth=1)
    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10
    plt.tick_params(axis='x',labelsize=10)
    plt.tick_params(axis='y',labelsize=10)

#    ax1.legend([title])
#    plt.xlim(0, nz-1)
    plt.xlabel(xlab, fontsize=10)
    plt.ylabel(ylab, fontsize=10)
    plt.xscale("log")
    plt.yscale("log")
#    plt.xlabel('log(Ra)', fontsize=10)
#    plt.ylabel('log(Nu)', fontsize=10)
    plt.legend(fontsize=10)

    full_filename = dir_path + filename + '.png'
    try: os.remove(full_filename)
    except OSError: pass
    plt.savefig(full_filename)
    plt.close()
    print("Saved file = ",full_filename)


#================================================================================
# Main
#================================================================================

HELP  = "no"
VALID = False
if len(argv) == 2:
    script, arg2 = argv
    if arg2 == "-h":
        HELP = arg2
        VALID = True
elif len(argv) == 3:
    script, PROJECT, RUN = argv
    VALID = True
else:
    f_project = open("PROJECT")
    f_run     = open("RUN")
    PROJECT = np.loadtxt(f_project)
    RUN     = np.loadtxt(f_run)

if not VALID or HELP == "-h":

#   Self documentation option

    print("============================================================================",file=sys.stderr,flush=True)
    if not VALID:
        print("ERROR: Incorrect usage of plot_bmk.py",file=sys.stderr,flush=True)
    else:
        print("plot_bmk.py documentation",file=sys.stderr,flush=True)
    print("============================================================================",file=sys.stderr,flush=True)

    print("Usage: python3 plot_bmk.py arguments (one to two arguments as below)",file=sys.stderr,flush=True)
    print(file=sys.stderr,flush=True)

    print("       python3 plot_bmk.py -h          (print documentation)",file=sys.stderr,flush=True)
    print("       python3 plot_bmk.py PROJECT RUN ",file=sys.stderr,flush=True)
    print("----------------------------------------------------------------------------",file=sys.stderr,flush=True)
    print("Note: PROJECT and RUN are the project & run names",file=sys.stderr,flush=True)
    print("      If unspecified, as in PROJECT and RUN files",file=sys.stderr,flush=True)


print(flush=True)
print("===============================================================",flush=True)
print("Making benchmark plot files for PROJECT = %s, RUN = %s" %(PROJECT,RUN),flush=True)
print("===============================================================",flush=True)
print(flush=True)
time.sleep(0.5)

# Set directory names

in_dir       = "./projects/" + PROJECT + "/"
plot_bmk_dir = "./projects/" + PROJECT +"/" + "/plots_bmk"

# Set the strings for the input and output file names

in_file    = in_dir   + RUN + ".in"
bmk_file   = in_dir   + RUN + ".bmk"

#==============================================================
# Check the directories exist and if necessary make them

#   Exit if the project input directory or input run file does not exist

if not os.path.exists(in_file):
    if rank == 0:
        print("ERROR: input file %s does not exist!" %(in_file),file=sys.stderr,flush=True)
    exit()

#   Exit if the project input directory and input run file does not exist

if not os.path.exists(bmk_file):
    if rank == 0:
        print("ERROR: benchmark file %s does not exist!" %(bmk_file),file=sys.stderr,flush=True)
    exit()

#   Make any required output directories

if rank == 0:
    if not os.path.exists(plot_bmk_dir):  os.mkdir(plot_bmk_dir)

# Get the parameters from the input file

import get_parameters
# Get the parameters for the run including plot choices and scales

get_parameters(in_file)
from get_parameters import *

if ny == 1:
    D = 2
else:
    D = 3

#-------------------------------------------------------------------------------
# Initialize for the desired geometry
#-------------------------------------------------------------------------------

# Initialize for the desired geometry

if geometry == "annular" or (geometry == "spherical" and D == 2):
    L_l,nz,R_p,L_p,radii,thetas,nu_factor = get_annular_geometry(R_p,L_p,theta1,theta2,nx,rank=0)
    if rank == 0:
        plot_annular_lattice(nz,radii,nx,theta1,theta2,annulus_interval)
elif geometry == "spherical" and D == 3:
    if rank == 0:
        print()
        print("ERROR: Spherical 3D geometry not yet implemented in TLBM_fast, exiting!",flush=True)
        print()
    exit()
else:
    R_p,L_p = get_R_and_L(R_p,L_p) # Get R_p and L_p from descriptive names if needed
    L_l     = nz-1

dx_p     = L_p/L_l       # Lattice spacing in physical units

#-------------------------------------------------------------------------------
# Calculate the lattice viscosity, thermal diffusivity, gravity and
# physical time step, and update Ra_p and Pr if unspecified
#-------------------------------------------------------------------------------

time_per_timestep, timescale_time_per_timestep, error = timevalue("time_per_timestep",time_per_timestep)
time_per_timestep = timevalue_to_secs(time_per_timestep,timescale_time_per_timestep) # Convert to secs if timescale != unspecified

Ra_p, Pr_p, eta_p, nu_p, ka_p, Ra_l, Pr_l, g_l, nu_l, ka_l, dt_p, dt_dt0 = get_lattice_viscosity(state, power_correction,Ra,Pr,g_p,beta_p,dT_p,L_p,eta_p,rho_p,ka_p,g_l,L_l,nu_l,nu_l_min,time_per_timestep,gamma_T,rheology,rank)

#-------------------------------------------------------------------------------
# Get nt, and the plot and data time increments
#-------------------------------------------------------------------------------

error = False
nt           , timescale_nt           , error = timevalue("nt"       ,nt       ,error)
Inc_data     , timescale_Inc_data     , error = timevalue("Inc_data" ,Inc_data ,error)

nt            = timevalue_to_timesteps(nt      ,timescale_nt       ,dt_p)
Inc_data      = timevalue_to_timesteps(Inc_data,timescale_Inc_data ,dt_p)

run_title   = str2str(run_title)
description = str2str(description)

rank = 0

#==============================================================
# Check the directories exist and if necessary make them

in_dir_path   = "./projects/bench2/"
dir_path      = "./bench_plots/"

vrms_max = 7.0
option = 8
if option == 1:
    labs = [ r'$Ra=10^4$', r'$Ra=10^5$', r'$Ra=10^6$', r'$Ra=10^7$'] 
    runs = ["Ra4Pr3", "Ra5Pr3", "Ra6Pr3", "Ra7Pr3"]
    BC_base = "ufs_dfx"
    title         = "BC: free slip (upper), fixed (lower)"
elif option == 2:
    labs = [ r'$Ra=10^4$', r'$Ra=10^5$', r'$Ra=10^6$', r'$Ra=10^7$', r'$Ra=10^8$', r'$Ra=10^9$']
    runs = ["Ra4Pr3_fs", "Ra5Pr3_fs", "Ra6Pr3_fs", "Ra7Pr3_fs", "Ra8Pr3_fs", "Ra9Pr3_fs"]
    BC_base       = "ufs_dfs"
    title         = "BC: free slip (upper & lower)"
elif option == 3:
    labs = [ r'$Ra=10^4$', r'$Ra=10^5$', r'$Ra=10^6$']
    cols = [ "blue"      , "green"     , "red"       ]
    runs = ["Ra4Pr1b", "Ra5Pr1b"       ,"Ra6Pr1b"    ]
    BC_base       = "pd_fx"
    Pr = 1
    title         = "BC: periodic (sides), fixed (upper/lower)"
    vrms = "physical"
    labs = [ r'$Ra=10^4,\ Pr=1$', r'$Ra=10^4,\ Pr=10$', r'$Ra=10^4,\ Pr=100$',  r'$Ra=10^4,\ Pr=1000$']
    cols = [ "blue"  , "green"   , "orange"  , "red"]
    runs = ["Ra4Pr1c", "Ra4Pr10c","Ra4Pr100c", "Ra4Pr1000c"   ]

    labs = [ r'$Ra=10^4,\ Pr=1000$', r'$Ra=10^5,\ Pr=1000$', r'$Ra=10^6,\ Pr=1000$', r'$Ra=10^7,\ Pr=1000$']
    runs = ["Ra4Pr1000f","Ra5Pr1000f","Ra6Pr1000f","Ra7Pr1000f"]

    runs = runs[:3]
elif option == 4:
    cols = [ "blue"  , "aqua", "orange"  , "red"]
    BC_base       = "pd_fx"
    title   = "BC: periodic (sides), fixed (upper/lower)"
    labs = [ r'$Ra=10^4,\ Pr=1000$', r'$Ra=10^5,\ Pr=1000$', r'$Ra=10^6,\ Pr=1000$',r'$Ra=10^7,\ Pr=1000$']
    runs = ["Ra4Pr1000","Ra5Pr1000","Ra6Pr1000","Ra7Pr1000"]
elif option == 5:
    in_dir_path   = "./projects/bench3/"
    cols = [ "blue"  , "aqua", "orange"  , "red"]
    BC_base       = "Pr1_fx"
    title   = "BC: fixed"
    labs = [ r'$Ra=10^4,\ Pr=1$', r'$Ra=10^5,\ Pr=1$', r'$Ra=10^6,\ Pr=1$',r'$Ra=10^7,\ Pr=1$']
    runs = ["Ra4Pr1fx","Ra5Pr1fx","Ra6Pr1fx","Ra7Pr1fx"]
elif option == 6:
    in_dir_path   = "./projects/bench_annular/"
    cols = [ "blue"  , "aqua", "orange"  , "red"]
    BC_base       = "annular_fx"
    title   = "BC: fixed"
    labs = [ r'$Ra=10^4$' + ", fixed", r'$Ra=10^5$' + ", fixed", r'$Ra=10^6$' + ", fixed",r'$Ra=10^7' + ",fixed"]
    runs = ["Ra4fx","Ra5fx","Ra6fx","Ra7fx"]
elif option == 7:
    in_dir_path   = "./projects/bench_annular/"
    cols = [ "blue"  , "aqua", "orange"  , "red"]
    BC_base       = "annular_fr"
    title   = "BC: free"
    labs = [ r'$Ra=10^4$' + ", free", r'$Ra=10^5$' + ", free", r'$Ra=10^6$' + ", free",r'$Ra=10^7' + ",free"]
    runs = ["Ra4fr","Ra5fr","Ra6fr","Ra7fr"]
elif option == 8:
    in_dir_path   = "./projects/bench_annular/"
    cols = [ "blue"  , "aqua", "orange"  , "red"]
    BC_base       = "annular_fr_Pr4"
    title   = "BC: free"
    labs = [ r'$Ra=10^4,\ Pr=10^4$' + ", free", r'$Ra=10^5,\ Pr=10^4$' + ", free", r'$Ra=10^6,\ Pr=10^4$' + ", free",r'$Ra=10^7,\ Pr=10^4' + ",free"]
    runs = ["Ra4Pr4fr","Ra5Pr4fr","Ra6Pr4fr","Ra7Pr4fr"]

R_p   = 6.378e6
L_p   = 2.9e6
factor_u = 1
factor_d = R_p/L_p
L_l = 511
dx_p = 2.9e6/L_l
sec_y = 365.25*24*3600
cm_y_2_m_s = 1/(100*sec_y)

Ra_Nu_base = "Ra_vs_Nu_" + BC_base

# Declare arrays as a function of Ra

n_ra = len(runs)

Ra_s      = np.zeros(n_ra)
Pr_s      = np.zeros(n_ra)
u_l_s     = np.zeros(n_ra)
vrms_l_s  = np.zeros(n_ra)
Re_l_s    = np.zeros(n_ra)
Re_vrms_l_s = np.zeros(n_ra)
Nu_pm_s   = np.zeros(n_ra)
Nu_gm_s   = np.zeros(n_ra)

# bmk file specifications

i_it     = 0
i_dt     = 1
i_dt_p   = 2
i_Ra     = 3
i_Pr     = 4
i_nu_l   = 5
i_u_l    = 6
i_vrms_l = 7
i_vrms_p = 8
i_vrms_cm_py = 9
i_Re_l   = 10
i_Re_vrms_l = 11
i_Re_p = 12
i_Nu_pm  = 13
i_Nu_gm  = 14

# Loop read the bmk files for each Ra and generate each "what" plot vs time

whats = ["Re_l","Nu_pm"]

for iwhat in range(len(whats)):
    what = whats[iwhat]

    fig = plt.figure(1,figsize=(4.5,3.5))
    plt.subplots_adjust(left=0.15, bottom=0.14, top=0.92, right=0.97)
    plt.title(title,fontsize=12)

    for i in range(n_ra):

        base_filename = runs[i]

        in_filename   = in_dir_path + base_filename + ".bmk"

        f_bmk = open(in_filename)
        bmk   = np.loadtxt(f_bmk)
        nt    = bmk.shape[0]

        itime     = np.zeros(nt)
        dt        = np.zeros(nt)
        dt_p      = np.zeros(nt)
        log_Ra    = np.zeros(nt)
        Ra        = np.zeros(nt)
        Pr        = np.zeros(nt)
        u_l       = np.zeros(nt)
        nu_l      = np.zeros(nt)
        vrms_l    = np.zeros(nt)
        Re_l      = np.zeros(nt)
        Re_vrms_l = np.zeros(nt)
        Nu_pm     = np.zeros(nt)
        Nu_gm     = np.zeros(nt)

        for it in range(nt):
            itime    [it] = bmk[it,i_it]
            dt_p     [it] = bmk[it,i_dt_p]
            log_Ra   [it] = bmk[it,i_Ra]
            Ra       [it] = 10**log_Ra[it]
            Pr       [it] = bmk[it,i_Pr]
            dt       [it] = bmk[it,i_dt]
            nu_l     [it] = bmk[it,i_nu_l]
            u_l      [it] = bmk[it,i_u_l]
            vrms_l   [it] = bmk[it,i_vrms_l]
            Re_l     [it] = bmk[it,i_Re_l]
            Re_vrms_l[it] = bmk[it,i_Re_vrms_l]
            Nu_pm    [it] = bmk[it,i_Nu_pm]
            Nu_gm    [it] = bmk[it,i_Nu_gm]

        it0 = 3*nt//4
        nt_ave = nt-it0-1
        dt_p        = dt_p[nt-1]
        Ra_s    [i] = Ra    [nt-1]
        Pr_s    [i] = Pr    [nt-1]
        u_l_s   [i] = u_l   [nt-1]/dt[nt-1]
        vrms_l_s[i] = vrms_l[nt-1]/dt[nt-1]
        Re_l_s  [i]    = np.sum(Re_l[it0:nt])     /nt_ave
        Re_vrms_l_s[i] = np.sum(Re_vrms_l[it0:nt])/nt_ave
        Nu_pm_s [i] = np.sum(Nu_pm[it0:nt])/nt_ave
        Nu_gm_s [i] = np.sum(Nu_gm[it0:nt])/nt_ave

        if   what == "vrms_l":
            f = vrms_l
        elif what == "u_l":
            f = u_l
        elif what == "Re_l":
            f = Re_l
        elif what == "Re_vrms_l":
            f = Re_vrms_l
        elif what == "Nu_pm":
            f = Nu_pm
        elif what == "Nu_gm":
            f = Nu_gm
        elif what == "dt":
            f = dt

#       Plot what vs time

        nt_out = nt
        plt.plot(itime[:nt_out]*dt_p/sec_my,f[:nt_out],"-",label=labs[i],color=cols[i],linewidth=1.5)

    plt.rcParams['xtick.labelsize'] = 10
    plt.rcParams['ytick.labelsize'] = 10
    plt.tick_params(axis='x',labelsize=10)
    plt.tick_params(axis='y',labelsize=10)
    if what == "Nu_pm" or what == "Nu_gm" or what == "Re_l" or what == "Re_vrms_l":
        plt.yscale("log")

#    plt.ylim(0)
    if False:
        if what == "Re_l":
            plt.ylim(0,5)
        if what == "Re_vrms_l":
            plt.ylim(0,2)
    plt.xlabel(r'My', fontsize=10)
    plt.ylabel(what, fontsize=10)
    plt.legend(fontsize=10)

    filename      = dir_path + what + "_" + BC_base
    full_filename = filename + '.png'
    try: os.remove(full_filename)
    except OSError: pass
    plt.savefig(full_filename)
    plt.close()
    print("Saved file = ",full_filename)

    Ra_s    [i] = Ra    [nt-1]
    Pr_s    [i] = Pr    [nt-1]
    u_l_s   [i] = u_l   [nt-1]
    vrms_l_s[i] = vrms_l[nt-1]
    Re_l_s  [i] = Re_l  [nt-1]
    Re_vrms_l_s[i] = Re_vrms_l[nt-1]
    Nu_pm_s [i] = Nu_pm [nt-1]
    Nu_gm_s [i] = Nu_gm [nt-1]

print("------------------------")

print("Ra_s     = ",Ra_s)
print("Pr_s     = ",Pr_s)
print("u_l_s    = ",u_l_s)
print("vrms_l_s = ",vrms_l_s)
print("Re_l_s   = ",Re_l_s)
print("Re_vrms_l_s = ",Re_vrms_l_s)
print("Nu_pm_s  = ",Nu_pm_s)
print("Nu_gm_s  = ",Nu_gm_s)

#--------------------------------------------------------------------------------

# Loop over each "what" to generate log(what) vs log(Ra)

whats = ["Re_l","Re_vrms_l", "Nu_pm"]

for iwhat in range(len(whats)):
    what = whats[iwhat]
    f = np.zeros(n_ra)
    for i in range(n_ra):
        if   what == "Re_l":
            f[i] = Re_l_s[i]
        elif what == "Re_vrms_l":
            f[i] = Re_vrms_l_s[i]
        elif what == "Nu_pm":
            f[i] = Nu_pm_s[i]
        elif what == "Nu_gm":
            f[i] = Nu_gm_s[i]

    c1,c0 = np.polyfit(np.log10(Ra_s), np.log10(f), 1)

    x_data = np.zeros(n_ra)
    y_data = np.zeros(n_ra)
    err    = np.zeros(n_ra)

    e2 = 0.0
    x_data = np.log10(Ra_s)
    y_data = c0 + c1*x_data
    err = y_data - np.log10(f)
    e2  = np.sum(err**2)
    stddev = np.sqrt(e2/n_ra)
    gamma = c1
    print("gamma = ",gamma,", Log_Ra vs Log_" + what + ", error = ",stddev)
    filename  = "log_Ra_vs_log_" + what + "_" + BC_base 
    plot_graph_line(title,filename,f,Ra_s,10**x_data,10**y_data,r'$Ra$',what)
