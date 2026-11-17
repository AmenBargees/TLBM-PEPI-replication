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
    plt.savefig(full_filename, bbox_inches='tight')
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
    with open("PROJECT") as f: PROJECT = f.read().strip()
    with open("RUN")     as f: RUN     = f.read().strip()

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

rank = 0

if not os.path.exists(bmk_file):
    print("ERROR: benchmark file %s does not exist!" %(bmk_file),file=sys.stderr,flush=True)
    exit()

os.makedirs(plot_bmk_dir, exist_ok=True)

in_dir_path = "./projects/" + PROJECT + "/"
dir_path    = "./projects/" + PROJECT + "/" + RUN + "/plots_bmk/"
os.makedirs(dir_path, exist_ok=True)

sec_y = 365.25*24*3600
sec_my = 1e6 * sec_y

BC_base = "runs"
title   = "PROJECT: " + PROJECT + "  RUN: " + RUN

# Read runs from RUNS file (multiple) or fall back to single RUN
if os.path.exists("RUNS"):
    with open("RUNS") as f:
        runs = f.read().split()
else:
    runs = [RUN]

cols = ["blue", "aqua", "orange", "red", "green", "purple", "brown", "pink"][:len(runs)]
labs = runs

# Label map for y-axis
what_labels = {
    "Re_vrms":   r'$Re_{vrms}$',
    "Re_vrms_l": r'$Re_{vrms,l}$',
    "Re_l":      r'$Re_l$',
    "Nu_pm":     r'$Nu$',
    "Nu_gm":     r'$Nu^*$',
    "vrms_l":    r'$v_{rms,l}$',
    "u_l":       r'$u_l$',
    "dt":        r'$\Delta t$',
}

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

        if not os.path.exists(in_filename):
            print(f"Skipping {base_filename}: {in_filename} not found")
            continue
        f_bmk = open(in_filename)
        bmk   = np.loadtxt(f_bmk)
        if bmk.ndim == 1: bmk = bmk[np.newaxis, :]
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
    plt.xlabel(r'$t$ (My)', fontsize=10)
    plt.ylabel(what_labels.get(what, what), fontsize=10)
    plt.legend(fontsize=10)

    filename      = dir_path + what + "_" + BC_base
    full_filename = filename + '.png'
    try: os.remove(full_filename)
    except OSError: pass
    plt.savefig(full_filename, bbox_inches='tight')
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

    mask = (Ra_s > 0) & (f > 0)
    if mask.sum() < 2:
        print(f"Skipping Ra scaling fit for {what}: not enough data points")
        continue
    c1,c0 = np.polyfit(np.log10(Ra_s[mask]), np.log10(f[mask]), 1)

    n_ok   = mask.sum()
    x_data = np.log10(Ra_s[mask])
    y_data = c0 + c1*x_data
    err    = y_data - np.log10(f[mask])
    e2     = np.sum(err**2)
    stddev = np.sqrt(e2/n_ok)
    gamma  = c1
    print("gamma = ",gamma,", Log_Ra vs Log_" + what + ", error = ",stddev)
    filename  = "log_Ra_vs_log_" + what + "_" + BC_base
    plot_graph_line(title,filename,f[mask],Ra_s[mask],10**x_data,10**y_data,r'$Ra$',what_labels.get(what,what))
