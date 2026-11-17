#===============================================================================
# Thermal Lattice Boltzmann code for mantle convection in 2D or 3D allowing
# for D2Q9 or D3Q19 using either single or multiple relaxation times, 1D or 2D
# domain decomposition and the option to use optimized c code for the TLBM
# kernel routines. See the README file for usage, and example input files
# archean.in, mantle.in, mantle3D and hadean.in for documentation on input parameters,
# which respectively show test runs of the transition from stagnant lid to
# plate tectonics in the archean, whole mantle convection in an annular model,
# whole mantle convection in a 3D model with an aspect ration 7.2:1 = 1/4 of
# the mantle, and turbulence in the Hadean.
#-------------------------------------------------------------------------------
# Authorship:
#
#     Peter Mora and Gabriele Morra (prototype python code)
#     Peter Mora (optimized c & python code, domain decomposition, optimal getting of edges,
#     data structure & cache blocking of c code, overlapping, graphics, unit conversion,
#     temperature dependent viscosity and plastic rheology as a yield stress power law,
#     input parameters and example input files, and scripts)
#     Gabriele Morra (Blankenbach benchmarking and miscellaneous)
#     Peter Mora (multiple relaxation time option - being re-implemented)
#
release      = "VERSION"
release_date = "RELEASE"
version_date = "DATE"
#-------------------------------------------------------------------------------
# Version 1 release notes:
#-------------------------------------------------------------------------------
# Version 1 is the first clean TLBM_fast version and can be used to model geodynamics scenarios with a
# wide range of Ra and Pr numbers (use Pr > 1000 for the present Earth), in 2D or 3D cartesian and
# annular models. It allows easy specification of plastic rheology temperature dependence, depth and/or
# spatial dependence of yield stress, temperature, viscosity and other key properties to be initialized
# as needed. Search for the keyword USERS to see where such initializations can be made.
#
# Version 1.01: Same as Version 1.0 except has some typos fixed.
# Version 1.10  Annular vector graphics fixed, better commenting of functions, better named options
#               for vector graphics, a more logical name for random temperature fluctuations (dT_p_random),
#               an option to initialize sinusoidal fluctuations which generate a specified number of
#               convection, and simpler & logical velocity filenames: u_mag to u and u_vectors to vectors,
#               and more typos fixed. Also, if the calculated lattice viscosity is too low (eg. < nu_l_min)
#               and needs to be increased, this changes the physical time step so the output increments are
#               automatically adjusted such that the output plots remain at Inc_plots*time_per_timestep.
#
# Current work is aimed at adding temperature and pressure dependent viscosity, partial melting,
# validating the LBM solutions for the 2D annular geometry option, developing a 3D spherical geometry option,
# testing the limits of scalability on HPC CPU clusters, and comparing results to those obtained with existing
# codes such as for modelling the transition from stagnant lid convection to plate tectonics in the Archean.
#
# An multiple relaxation time option from earlier prototype codes will be remimplemented
# in a future version, to add to the current single relaxation time cod
#-------------------------------------------------------------------------------
# COPYRIGHT
#
# Strictly no distribution allowed outside of the collaboration between
# Peter Mora and Gabriele Morra without permission from both authors
#===============================================================================

from ctypes import *
from mpi4py import MPI
import os
from os.path import exists
import matplotlib
import numpy as np
import math
import matplotlib.pyplot as plt
import random
from numpy.linalg import inv
from numpy.random import randn, SeedSequence, default_rng
from mpl_toolkits.mplot3d import axes3d, Axes3D
import time
import timeit
import sys
from sys import argv
from datetime import datetime
import subprocess 
from scipy.interpolate import interp1d
#===============================================================================
# Python LBM kernels
#===============================================================================

#-------------------------------------------------------------------------------
# Python multiple relaxation time collision step
#-------------------------------------------------------------------------------

def collision_MRT_py():
    global rho, u, f, tau_f, Mi, Mi_S, S_, w, shifted_distribution

    m_eq    = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float)
    if D == 2:
        jx      = rho * u[:,:,X,:]
        jy      = rho * u[:,:,Y,:]
        jx2     = jx * jx
        jy2     = jy * jy
        j2      = jx2 + jy2
        threej2 = 3 * j2
        m_eq[:,:,0,:] = rho
        m_eq[:,:,1,:] = -2 * rho + threej2
        m_eq[:,:,2,:] = rho - threej2
        m_eq[:,:,3,:] =   jx
        m_eq[:,:,4,:] = - jx
        m_eq[:,:,5,:] =   jy
        m_eq[:,:,6,:] = - jy
        m_eq[:,:,7,:] =   jx2 - jy2
        m_eq[:,:,8,:] =   jx  *  jy
    elif D == 3:
        jx       = rho * u[:,:,X,:]
        jy       = rho * u[:,:,Y,:]
        jz       = rho * u[:,:,Z,:]
        jx2      = jx * jx
        jy2      = jy * jy
        jz2      = jz * jz
        j2       = jx2 + jy2 + jz2
        m_eq[:,:, 0,:] = rho
        m_eq[:,:, 1,:] = -11 * rho + 19  * j2
        m_eq[:,:, 2,:] = 3 * rho   - 5.5 * j2
        m_eq[:,:, 3,:] = jx
        m_eq[:,:, 4,:] = -0.66666666666667 * jx
        m_eq[:,:, 5,:] = jy
        m_eq[:,:, 6,:] = -0.66666666666667 * jy
        m_eq[:,:, 7,:] = jz
        m_eq[:,:, 8,:] = -0.66666666666667 * jz
        m_eq[:,:, 9,:] = 2 * jx2 - jy2 - jz2
        m_eq[:,:,10,:] = -(0.5) * m_eq[:,:,9,:]
        m_eq[:,:,11,:] = jy2 - jz2
        m_eq[:,:,12,:] = -(0.5) * (jy2 - jz2)
        m_eq[:,:,13,:] = jx * jy
        m_eq[:,:,14,:] = jy * jz
        m_eq[:,:,15,:] = jx * jz
        m_eq[:,:,16,:] = 0
        m_eq[:,:,17,:] = 0
        m_eq[:,:,18,:] = 0

    m       = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float)
    for a in range(na):
        for b in range(na):
            if shifted_distribution:
                m[:,:,a,:] += M[a,b] * (f[:,:,b,:]+w[b])
            else:
                m[:,:,a,:] += M[a,b] *  f[:,:,b,:]

    Delta_m = m_eq - m
    omega = 1/tau_f         # This gets the desired viscosity
    Delta_f = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float)
    for a in range(na):
        for b in range(na):
            if D == 2:
                if b < 7:
                    Delta_f[:,:,a,:] += (Mi_S[a,b] * Delta_m[:,:,b,:])
                else:
                    Delta_f[:,:,a,:] += (Mi[a,b] * omega * Delta_m[:,:,b,:])
            elif D == 3:
                if b < 9 or b == 10 or b == 12 or b > 15:
                    Delta_f[:,:,a,:] += (Mi_S[a,b] * Delta_m[:,:,b,:])
                else:
                    Delta_f[:,:,a,:] += (Mi[a,b] * omega  * Delta_m[:,:,b,:])
    return Delta_f

#-------------------------------------------------------------------------------
#   Print the status if the file was successfully read
#-------------------------------------------------------------------------------

    if rank == 0:
        for iout in range(2):
            if   iout == 0: output_print = sys.stderr
            elif iout == 1: output_print = sys.stdout
            if nwarning > 0:
                print(" file read with %d warnings." %(nwarning),file=output_print,flush=True)
            else:
                print(" file read successfully!",file=output_print,flush=True)
            print("------------------------------------------------------------------",file=output_print,flush=True)

#===============================================================================
# Function to plot an annular lattice
#===============================================================================

def plot_annular_lattice(n_r,radii,n_theta,theta1,theta2,annulus_interval):

#   Calculate the plot resolution

    theta1 *= np.pi/180  ;theta2 *= np.pi/180

    n_r_new  = n_r    //annulus_interval 
    if n_r_new < 3:  annulus_interval = max(n_r//3,1)
    n_r      = n_r    //annulus_interval
    n_theta  = n_theta//annulus_interval

    theta     = np.linspace(theta1,theta2,n_theta)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)

#   Calculate r, x and y for the annular lattice plot

    r         = np.zeros(n_r)
    for i_r in range(n_r):
        r[i_r] = radii[i_r*annulus_interval]/1000       # Do the plot in km, not m

    x         = np.zeros((n_r,n_theta))
    y         = np.zeros((n_r,n_theta))

    for i_r in range(n_r):
        x[i_r,:] = r[i_r]*cos_theta
        y[i_r,:] = r[i_r]*sin_theta

    x_min = np.min(x) ; x_max = np.max(x)
    y_min = np.min(y) ; y_max = np.max(y)
    aspect = (x_max-x_min)/(y_max-y_min)

    plot_size0  = 2 + aspect/2

    if aspect < 1.5: plot_size0 *= 1.2

    space_left   = 0.57
    space_bottom = 0.39
    space_top    = 0.24
    space_right  = 0.1

    if   plot_size == "vsmall": plot_size0 *= 0.8
    elif plot_size == "small" : plot_size0 *= 1.0
    elif plot_size == "medium": plot_size0 *= 1.5
    elif plot_size == "large" : plot_size0 *= 2.0
    elif plot_size == "vlarge": plot_size0 *= 2.5
    else:                       plot_size0 = float(plot_size)

    xsize        = plot_size0
    ysize        = plot_size0/aspect

    xsize        = plot_size0        + space_left + space_right
    ysize        = plot_size0/aspect + space_top  + space_bottom

#   Space as a fraction of xsize and ysize for adjusting space by subplots_adjust

    space_left   /= xsize
    space_right  /= xsize
    space_bottom /= ysize 
    space_top    /= ysize 

    plot_size0  *= plot_size_factor

    plt.figure(figsize=(xsize,ysize))
    fig = plt.subplots_adjust(left=space_left,bottom=space_bottom,top=1-space_top,right=1-space_right)

    if plot_size == "vsmall":
        fontsize_title       = 7
        fontsize_axis_labels = 6
        fontsize_axis_tics   = 5
    elif plot_size == "small":
        fontsize_title       = 9
        fontsize_axis_labels = 7
        fontsize_axis_tics   = 6
    elif plot_size == "medium":
        fontsize_title       = 10
        fontsize_axis_labels = 9
        fontsize_axis_tics   = 8
    else:
        fontsize_title       = 12
        fontsize_axis_labels = 10
        fontsize_axis_tics   = 8

    plt.title("Annular lattice",fontsize=fontsize_title)
    for i_theta in range(n_theta):
        plt.plot(x[:,i_theta],y[:,i_theta],"k-",linewidth=0.5)
    for i_r in range(n_r):
        plt.plot(x[i_r,:],y[i_r,:],"k-",linewidth=0.5)

    plt.xlabel("x (km)",fontsize=fontsize_axis_labels)
    plt.ylabel("z (km)",fontsize=fontsize_axis_labels)
    plt.tick_params(axis='x',labelsize=fontsize_axis_tics)
    plt.tick_params(axis='y',labelsize=fontsize_axis_tics)

    full_filename = plot_dir + "annular_lattice.png"

    try: os.remove(full_filename)
    except OSError: pass
    plt.savefig(full_filename)
    plt.close()

    print("Saved annular lattice to file ",full_filename)
    print("==================================================================",flush=True)

#===============================================================================
# Function to plot a field profile (ie. mean function of z)
#===============================================================================


def bmk_diagnostics():
    """V2.96-style benchmark diagnostics (Peter Mora's new conventions):
    finite-difference Nusselt with cylindrical conductive normalization,
    dx^2-weighted volumetric mean temperature, and volume-weighted
    nondimensional Vrms with the per-row conformal velocity conversion
    u_p ~ u_l * dx(z). COLLECTIVE: call on every rank; returns a dict on
    rank 0 and None elsewhere."""
    u2_loc = np.einsum('yzdx,yzdx->yzx', u, u)
    T_loc  = eta/DR_2
    sum_u2 = np.zeros(nz, dtype=np.float64)
    sum_T  = np.zeros(nz, dtype=np.float64)
    cnt    = np.zeros(nz, dtype=np.float64)
    z0 = Rz*nzi
    for zl in range(mz[Rz]):
        zg = z0 + zl
        sl_u = u2_loc[yedge:ny_i-yedge if yedge else ny_i, zedge+zl, xedge:nx_i-xedge if xedge else nx_i]
        sl_T = T_loc [yedge:ny_i-yedge if yedge else ny_i, zedge+zl, xedge:nx_i-xedge if xedge else nx_i]
        sum_u2[zg] = np.sum(sl_u, dtype=np.float64)
        sum_T [zg] = np.sum(sl_T, dtype=np.float64)
        cnt   [zg] = sl_u.size
    SUM_U2 = np.zeros_like(sum_u2) ; SUM_T = np.zeros_like(sum_T) ; CNT = np.zeros_like(cnt)
    comm.Reduce(sum_u2, SUM_U2, root=0)
    comm.Reduce(sum_T , SUM_T , root=0)
    comm.Reduce(cnt   , CNT   , root=0)
    if rank != 0: return None
    CNT   = np.maximum(CNT, 1.0)
    Tz    = SUM_T/CNT                           # centered T profile, index 0 = surface
    u2z   = SUM_U2/CNT
    Tz_nd = (Tz + dT_p/2)/dT_p                  # nondimensional [0,1]
    if geometry == "annular":
        r_nd = radii/np.float64(L_p) ; dxp = (radii/radii[0]).astype(np.float64)
        r2   = float(r_nd[0]) ; r1 = float(r_nd[-1])
        q_cond = (1.0/r2)/np.log(r2/r1)
    else:
        dxp  = np.ones(nz) ; q_cond = 1.0 ; r2 = 1.0
        r_nd = 1.0 - np.arange(nz)/(nz - 1.0)   # pseudo radius: plain gradient
    w  = dxp**2 ; wt = w.copy() ; wt[0] *= 0.5 ; wt[-1] *= 0.5
    T_bar    = float((Tz_nd*wt).sum()/wt.sum())
    T_mean_C = float((Tz*wt).sum()/wt.sum()) + T0_p + dT_p/2
    dTdr = (Tz_nd[1:] - Tz_nd[:-1])/(r_nd[1:] - r_nd[:-1])
    r_fd = 0.5*(r_nd[1:] + r_nd[:-1])
    q_fd = -(1.0/q_cond)*(r_fd/r2)*dTdr
    iu = 0
    for i in range(nz//2):
        if q_fd[i+1] < q_fd[i]:
            iu = i ; break
    idn = nz - 3
    for idn in range(nz-3, nz//2, -1):
        if q_fd[idn] < q_fd[idn+1]: break
    idn += 1
    Nu_u = float(q_fd[iu]) ; Nu_d = float(q_fd[idn]) ; Nu = 0.5*(Nu_u + Nu_d)
    v2w      = u2z*dxp**2                       # conformal conversion (u_p ~ u_l dx)
    Vrms_lat = math.sqrt(max(float((v2w*wt).sum()/wt.sum()), 0.0))
    Vrms_nd  = Vrms_lat*(L_l/ka_l)/dt           # /dt: vary-dt rescaling of lattice u
    v_rms_cm = Vrms_lat*dx_p/(dt*dt_p)*u_cm_py
    return {"Vrms_nd": Vrms_nd, "v_rms_cm": v_rms_cm, "Nu": Nu, "Nu_u": Nu_u,
            "Nu_d": Nu_d, "T_bar": T_bar, "T_mean_C": T_mean_C}


def plot_field_profile(field,what,what_latex,c0=0,c1=1):

#   Get the field, and then field profile on core 0

    FIELD = np_float(get_field(field))

    _bd = bmk_diagnostics() if what == "T" else None

    if rank == 0:
        FIELD = c0 + c1 * FIELD
        FIELD_z = np.zeros(nz,dtype=np_float)
        for z in range(nz):
            FIELD_z[z] = np.sum(FIELD[:,z,:])/(nx*ny)
 
#       Calculate the thermal boundary layer width and print it and Nu

        nz_profile = nz
        zu_TBL = np.zeros
        if geometry == "annular":
            z    = np.zeros(nz,dtype=np_float)
            for i in range(nz):
                z[:] = R_p - radii[:]
        else:
            z  = dx_p*np.linspace(0,nz-1,nz,dtype=np_float)
        if (write_data or write_plots) and what == "T":
            i0     = 1 # Start point to calculate TBL width (use 1 as the boundary point is not part of the linear region)
            if geometry == "annular":
                nu_TBL = n_TBL
                nd_TBL = n_TBL
            else:
                nu_TBL = n_TBL
                nd_TBL = n_TBL
            zu_TBL      = np.zeros(nz,dtype=np_float)
            zu_TBL[:nz] = z      [:nz]
            Tu_z        = FIELD_z[:nz]
            zd_TBL      = np.zeros(nz,dtype=np_float)
            Td_z        = np.zeros(nz,dtype=np_float)
            for i in range(nz):
                zd_TBL[i] = L_p - z[nz-i-1]
                Td_z  [i] = FIELD_z[nz-i-1]

            c1_Tu,c0_Tu = np.polyfit(zu_TBL[i0:nu_TBL+i0], Tu_z[i0:nu_TBL+i0], 1)
            c1_Td,c0_Td = np.polyfit(zd_TBL[i0:nd_TBL+i0], Td_z[i0:nd_TBL+i0], 1)

            w1_TBL = (T0_p+dT_p/2-c0_Tu)/c1_Tu if c1_Tu != 0 else np.inf
            w2_TBL = (T0_p+dT_p/2-c0_Td)/c1_Td if c1_Td != 0 else np.inf
            w_TBL  = (w1_TBL+w2_TBL)/2
            for i in range(nz):
                if zu_TBL[i] > w1_TBL: break
            nzu_profile = i+1
            for i in range(nz):
                if zd_TBL[i] > w2_TBL: break
            nzd_profile = i+1

            FIELDu_z = FIELD_z[:nz]
            FIELDd_z = np.zeros(nz,dtype=np_float)
            for i in range(nz): FIELDd_z[i] = FIELD_z[nz-i-1]
            Tu_pred  = c0_Tu + c1_Tu * zu_TBL
            Td_pred  = c0_Td + c1_Td * zd_TBL

#           Nu by Peter's   method of fitting and averaging upper and lower TBL widths

            Nu_pm = L_p/w_TBL

#           Nu by Gabriele's method of dinite difference one up from the lower surface

            dT_dz = FIELD[:,1:nz,:] - FIELD[:,0:nz-1,:]
            Nu_gm = (1/dT_p)*np.sum(dT_dz[0,nz-3,:])    # Original as per paper but allowing for nz != nx
#            dT_dz = FIELD_z[1:nz]-FIELD_z[0:nz-1]
#            dT_dz = dT_dz[nz-3]
#            Nu_gm = (1/dT_p)*dT_dz*(nz-1)

            Nu    = Nu_pm

#           V_rms physical and renormalized time

            v_rms_p       = v_rms_l*dx_p/dt_p
            v_rms_cm_py   = v_rms_p*u_cm_py
            Re_l          = u_mag_max*L_l/nu_l
            Re_rms_l      = v_rms_l*L_l/nu_l
            Re_p          = u_mag_max*L_p/(ka_p*Pr)
            Re_rms_p      = v_rms_p*L_p/(ka_p*Pr)
            u_l           = u_mag_max
            if _bd:
                print("BMK t = %s  Vrms_nd=%9.4g  v_rms=%8.3f cm/yr  Nu=%7.3f (u %6.3f, d %6.3f)  T_bar=%7.4f  T_mean=%7.1f C" % (
                    timesteps_to_physical_time_text(t,dt_p,output_units), _bd["Vrms_nd"], _bd["v_rms_cm"], _bd["Nu"], _bd["Nu_u"], _bd["Nu_d"], _bd["T_bar"], _bd["T_mean_C"]), flush=True)

            if output_units == "physical":
                time_text = timesteps_to_physical_time_text(t,dt_p,output_units)
                print("Time = %d (%d computed) = %s, " %(it_mach,t_computed,time_text),end="",flush=True)
            else:
                print("Time = %d (%d computed), " %(it_mach,t_computed),end="",flush=True)

#           Write benchmark file with: log(Ra) L_p wu_TBL wd_TBL V_rms

            f_bmk = open(bmk_file,"a")       # Open  benchmark file
            f_bmk.writelines("%d %11.4f %5.3f %9.4g %8.5f %8.4f %8.4f %8.4f %10.5g %10.4g %8.4g\n" % (
                it, it*dt_p/sec_my, math.log(Ra_l,10.0), dt, _bd["T_bar"], _bd["Nu"], _bd["Nu_u"], _bd["Nu_d"],
                _bd["Vrms_nd"], _bd["v_rms_cm"], mach_max))
            f_bmk.close()                    # Close benchmark file

        if write_plots:

            if not latex:
                title      = what
            else:
                title      = what_latex

#           Set and plot the title

            if run_title == "unspecified": this_title = PROJECT + "/" + RUN
            else:                          this_title = run_title

            this_title = str2str(this_title)

            if t < 0:
                title = this_title + ": " + title
            elif output_units == "physical":
                time_title = timesteps_to_physical_time_text(t,dt_p,output_units)
                title = this_title + ": " + title + "_z(%s)" %(time_title)
            else:
                title = this_title + ": " + title + "_z(%d)" %(it_plots)

            if what == "T":
                title = title + ", Nu=%4.2f" %(Nu)

#           Set the filename and save the figure

            what_file = what + "_z"

            if t < 0: full_filename = plot_dir + what_file + ".png"
            else:     full_filename = plot_dir + what_file + "_%s.png" %(timetext(it_plots,timescale_Inc_plots,dt_p))

            print("Plotting profile for %s E [%8.3g,%8.3g] to %s" %(what,np.min(FIELD_z),np.max(FIELD_z),full_filename),flush=True)

            fig = plt.figure(figsize=(4.7,3.6))
            fig = plt.subplots_adjust(left=0.15,bottom=0.12,top=0.935,right=0.98)
            plt.title(title,fontsize=9)
            labelu = what + "_upper"
            labeld = what + "_lower"
            plt.plot(zu_TBL,FIELDu_z,"-o",color="blue",label=labelu,markersize=1,linewidth=0.75)
            plt.plot(zd_TBL,FIELDd_z,"-o",color="red" ,label=labeld,markersize=1,linewidth=0.75)
            if what == "T":
                nz_plot = min(nz_profile, int(w_TBL + 2) if np.isfinite(w_TBL) else nz_profile)
                plt.plot(zu_TBL[i0:nzu_profile],Tu_pred[i0:nzu_profile],"--",color="blue",label="Fit upper",linewidth=0.6)
                plt.plot(zd_TBL[i0:nzd_profile],Td_pred[i0:nzd_profile],"--",color="red" ,label="Fit lower",linewidth=0.6)
                plt.axhline(T0_p+dT_p/2,linestyle="-",linewidth=0.5,color="black")
                plt.xlim(0,L_p/2)
                plt.ylim(0,T0_p+dT_p*1.2)
            plt.legend(loc="upper right",fontsize=9)
            plt.xlabel("|z-z_boundary|",fontsize=9)
            plt.ylabel(what,fontsize=9)
            plt.tick_params(axis='x',labelsize=8)
            plt.tick_params(axis='y',labelsize=8)

            try: os.remove(full_filename)
            except OSError: pass
            plt.savefig(full_filename)
            plt.close()

#===============================================================================
# Function to write a data file and plot a field
#===============================================================================

def write_and_plot_field(field,what,what_latex,cbar,cbar_latex,save=True,plot=True, \
    clip=1,gpow=1,zero_mean=False,plot_image=True,plot_vectors=False,vmin=0,vmax=0,c0=0,c1=1):

    global coordinate

    if save and not globals().get("_save_fields_now", True): save = False

    if (write_data and save) or (write_plots and plot):

        FIELD = get_field(field)

        time_text_print = timesteps_to_physical_time_text(it,dt_p,output_units)
        time_text_file  = timetext(it,timescale_Inc_data,dt_p)

        filename = data_dir + what + "_%s.npy" %(time_text_file)

        if write_data and save and rank == 0:
            FIELD_min = np.min(FIELD)
            FIELD_max = np.max(FIELD)
            print("Saving %s E [%9.4g,%9.4g]" %(what,FIELD_min,FIELD_max),end="",flush=True)
            np.save(filename,FIELD)
            if print_filename:
                print(", %s saved!" %(filename),flush=True)
            else:
                print()

        if write_plots and plot and plot_field and rank == 0:

            # Always use the LaTeX-style symbol for the plot title (e.g.
            # ρ_T) — matplotlib's built-in MathText renders $...$ markup
            # without needing system LaTeX. When latex=False (the
            # default, no usetex), MathText does NOT understand the
            # \displaystyle command, so we strip it; the math symbol
            # itself still renders fine.
            if not latex:
                title      = what_latex.replace(r'\displaystyle','').replace(r'\displaystyle{}','').strip()
                cbar_label = cbar
            else:
                title      = what_latex
                cbar_label = cbar_latex

            for i in range(len(coordinates)):

                coordinate = int(coordinates[i])
                FIELDi     = get_plane(FIELD,plane,coordinate)

                if output_units != "lattice":
                    FIELDi = convert_units(FIELDi,FIELDi.dtype,c0,c1)
                else:
                    if not latex: cbar_label = "Lattice units"
                    else:         cbar_label = r'Lattice units'

                plot_field(title,what,FIELDi,clip=clip,gpow=gpow,zero_mean=zero_mean,cbar_label=cbar_label,plot_image=plot_image,plot_vectors=plot_vectors,vmin=vmin,vmax=vmax)

#===============================================================================
# Function to plot a field
#===============================================================================
# There is an option to plot an image, and to plot_vectors (velocity vectors), so
# the velocity vectors can overlay an image. If an image is plotted, then the
# color of the vectors is as specified by vector_color in the input file.
# If plot_image = False, then the velocity vectors are black
#---------------------------------------------------------------------------------------

def plot_field(title,what,field,clip,gpow,zero_mean=False,cbar_label="",plot_image=True,plot_vectors=False,vmin=0,vmax=0):

    global vector_color, full_filename, coordinate

    field     = np.float64(field)
    field_min = np.min(field)
    field_max = np.max(field)

    print("Plotting ",title," E [%9.4g" %(field_min),",%9.4g" %(field_max),"]",end="",flush=True)
    if D == 3: print(", plane = %s, coordinate = %d" %(plane,coordinate),end="",flush=True)

#   If zero_mean = False, set the color range based on vmin and vmax if specified and
#   otherwise on the max and min values in the field. If zero_mean = True,
#   then the field is assumed to be a zero mean field like velocity so the range
#   will be from -max(|field|) to max(|field|) so field=0 will have the same color

    if vmin == 0 and vmax == 0:
        if zero_mean:
            field_abs_max = np.max(field_abs)
            vmin          = -field_abs_max
            vmax          =  field_abs_max
        else:
            vmin          = field_min
            vmax          = field_max

    field0 = (vmin+vmax)/2

#   Clip the field

    vmin = field0 - clip*(field0-vmin)
    vmax = field0 + clip*(vmax-field0)
    
#   Do gpow nonlinear scaling

    tiny                 = np.finfo(field.dtype).tiny
    dfield               = field - field0
    field_abs            = np.abs(dfield)
    max_dfield           = np.max(dfield)
    field_abs_normalized = field_abs/(np.max(field_abs) + 2*tiny)
    field                = field0 + max_dfield*np.sign(dfield)*np.power(field_abs_normalized,gpow)

#   Set the aspect ratio of the plot

    ny = field.shape[0]
    nx = field.shape[1]
    if geometry == "annular":
        dtheta = thetas[1] - thetas[0]
        theta = np.linspace(np.pi*theta1/180,np.pi*theta2/180+dtheta,nx)
        theta, radius = np.meshgrid(theta,radii)
        cos_theta = np.cos(theta)
        sin_theta = np.sin(theta)
        x_lim = radius*cos_theta
        y_lim = radius*sin_theta
        x_min = np.min(x_lim) ; x_max = np.max(x_lim)
        y_min = 0             ; y_max = np.max(y_lim)
        aspect = (x_max-x_min)/(y_max-y_min)
        if np.max(thetas)-np.min(thetas) <= 1.01*np.pi: 
             aspect = 1.5
        else:
             aspect = 1
        aspect0  = aspect
        aspect   = 1
    else:
        aspect = nx/ny

#   Calculate the plot size and set space required for the title,
#   axis labels and colorbar if it is to be plotted

    if geometry == "cartesian" or (geometry == "annular" and plot_image):
        if plot_size == "vsmall":
            fontsize_title       = 7
            fontsize_axis_labels = 6
            fontsize_axis_tics   = 5
        elif plot_size == "small":
            fontsize_title       = 9
            fontsize_axis_labels = 7
            fontsize_axis_tics   = 6
        elif plot_size == "medium":
            fontsize_title       = 10
            fontsize_axis_labels = 9
            fontsize_axis_tics   = 8
        else:
            fontsize_title       = 12
            fontsize_axis_labels = 10
            fontsize_axis_tics   = 8
    elif geometry == "annular":
        if plot_size == "vsmall":
            fontsize_title       = 6
            fontsize_axis_labels = 5
            fontsize_axis_tics   = 4
        elif plot_size == "small":
            fontsize_title       = 8
            fontsize_axis_labels = 6
            fontsize_axis_tics   = 5
        elif plot_size == "medium":
            fontsize_title       = 9
            fontsize_axis_labels = 7
            fontsize_axis_tics   = 6
        else:
            fontsize_title       = 12
            fontsize_axis_labels = 10
            fontsize_axis_tics   = 8

    space_left   = 0.55
    space_bottom = 0.39
    space_top    = 0.24

    if plot_colorbar and plot_image: space_right  = 0.6
    elif plot_vectors and geometry == "annular":
                                     space_top    = 0.25
                                     space_bottom = 0.1
                                     space_left   = 0.1
                                     space_right  = 0.1
    elif plot_vectors:               space_right  = 0.35
    else:                            space_right  = 0.1

    if geometry == "annular":
        if plot_colorbar and plot_image:
            space_bottom = 0.0
            space_left   = 0.1
            space_top   -= 0.1
            if   plot_size == "small":
                space_right -= 0.1
            elif plot_size == "vsmall":
                space_right -= 0.15
        elif plot_image:
            space_right -= 0.1
            space_bottom = 0.05
            space_left   = 0.05
            if plot_size == "vsmall": space_top   -= 0.05

    plot_size0 = 2 + aspect/2
    xinch0     = plot_size0

    if   plot_size == "vsmall": plot_size0 *= 0.8
    elif plot_size == "small" : plot_size0 *= 1.0
    elif plot_size == "medium": plot_size0 *= 1.5
    elif plot_size == "large" : plot_size0 *= 2.0
    elif plot_size == "vlarge": plot_size0 *= 2.5
    else:                       plot_size0 = float(plot_size)

    plot_size0  *= plot_size_factor

    xinch0 = plot_size0
    xinch  = plot_size0
    yinch  = plot_size0/aspect

    xsize = xinch + space_left + space_right
    ysize = yinch + space_top  + space_bottom

    if geometry == "annular":
        if np.max(thetas-np.min(thetas)) <= 1.01*np.pi:
            aspect_cbar       = 40*ysize/xsize      # Nice aspect ratio for the colorbar
            fraction_colorbar = 0.015               # Nice colorbar size fraction for annular geometry
        else:
            aspect_cbar       = 50*ysize/xsize      # Nice aspect ratio for the colorbar
            fraction_colorbar = 0.02                # Nice colorbar size fraction for annular geometry
    else:
        aspect_cbar       = 50*ysize/xsize          # Nice aspect ratio for the colorbar
        if aspect <= 1.5: fraction_colorbar = 0.022 # Nice colorbar size fraction
        else:             fraction_colorbar = 0.03  # Nice colorbar size fraction

    pad               = 0.2/xsize                   # Nice padding for colorbar

#   Allow extra space if the description is in the frame

    if len(description) > 0 and description_location == "frame":
        if   geometry == "cartesian":
            if aspect <= 4.0:
                space_bottom += 0.15
            else:
                space_bottom += 0.05
        elif geometry == "annular":
             if aspect0 < 1.25:
                 space_bottom += 0.2
                 space_top    += 0.1

#   Space as a fraction of xsize and ysize for adjusting space by subplots_adjust

    space_left   /= xsize
    space_right  /= xsize
    space_bottom /= ysize 
    space_top    /= ysize 

#   Set the figure size and plot 

    fig = plt.figure(figsize=(xsize,ysize))

    fig.subplots_adjust(left=space_left,bottom=space_bottom,top=1-space_top,right=1-space_right)
    if geometry != "annular":
        fig = plt.subplots_adjust(left=space_left,bottom=space_bottom,top=1-space_top,right=1-space_right)
        x_min = 0
        z_max = 0
        if output_units == "lattice":
            x_max  = nx-1
            z_min  = ny-1
            factor = 1
        elif L_p < 1000000:            # Small models like upper mantle => metres
            x_max  = dx_p*(nx-1)/1.e3
            z_min  = dx_p*(ny-1)/1.e3
            factor = dx_p/1.e3
        else:                          # Large models like full mantle => km
            x_max  = dx_p*(nx-1)/1.e6
            z_min  = dx_p*(ny-1)/1.e6
            factor = dx_p/1.e6

#       Plot the image if required

        if plot_image:
            image = plt.imshow(field,extent=[x_min,x_max,z_min,z_max],aspect=1,cmap=colormap,vmin=vmin,vmax=vmax)

    else:
        ax = fig.add_axes([space_left,space_bottom,1-space_left-space_right,1-space_top-space_bottom],polar=True)
        if plot_image:
            plt.xlim(theta1*np.pi/180,theta2*np.pi/180)
            plt.ylim(0,np.max(radii))
            image = ax.pcolormesh(theta,radii,field,edgecolors='face',cmap=colormap,vmin=vmin,vmax=vmax)
        ax.set_xticks([])       # No theta  ticks
        ax.set_yticklabels([])  # No radial tick  values

#   Overlay velocity vectors if required

    if plot_vectors:
        interval = vector_interval

        if plot_image == False: col_vectors = "k"
        else:                   col_vectors = vector_color
        if geometry == "cartesian":
            Y, X = np.mgrid[0:field.shape[0],0:field.shape[1]]*factor
            if   plane == "xz":
                U1 = U_x ; U2 = U_z
            elif plane == "xy":
                U1 = U_x ; U2 = U_y
            elif plane == "yz":
                U1 = U_y ; U2 = U_z

            u_mag = np.max(np.sqrt(U1**2+U2**2))
            uscale   = 4*xinch0/(u_mag+4*xinch0*tiny)
            length   = (nx-1)/interval

            print(", overlaying velocity vectors",end="",flush=True)

            plt.quiver(X[::interval,::interval],Y[::interval,::interval],U1[::interval,::interval]*uscale,-U2[::interval,0::interval]*uscale,color=col_vectors,width=0.002*vector_width,scale=length/vector_length,scale_units='inches')
        elif geometry == "annular":
            theta, radius = np.meshgrid(thetas,radii)
            U1 = (-U_z*np.cos(theta) - U_x*np.sin(theta))
            U2 = (-U_z*np.sin(theta) + U_x*np.cos(theta))

            u_mag = np.max(np.sqrt(U1**2+U2**2))
            uscale   = xinch0*L_p/R_p/(u_mag+2*L_p/R_p*xinch0*tiny)
            length   = (nz-1)/interval

            print(", overlaying velocity vectors",end="",flush=True)

            ax.quiver(theta[::interval,::interval],radius[::interval,::interval],U1[::interval,::interval]*uscale,U2[::interval,0::interval]*uscale,color=col_vectors,width=0.002*vector_width,scale=30/vector_length)

    if what == "u": print(flush=True)

#   Set and plot the title

    if run_title == "unspecified": this_title = RUN
    else:                          this_title = run_title

    this_title = str2str(this_title)

    if t < 0:
        title = this_title + ": " + title
    elif output_units == "physical":
        time_title = timesteps_to_physical_time_text(t,dt_p,output_units)
        title = this_title + ": " + title + "(%s)" %(time_title)
    else:
        title = this_title + ": " + title + "(%d)" %(it_plots)

    if D == 3:
        title = title + ", plane=" + plane
        if   plane == "xz": title = title + ", y=" + str(coordinate)
        elif plane == "xy": title = title + ", z=" + str(coordinate)
        elif plane == "yz": title = title + ", x=" + str(coordinate)

    plt.title(title,fontsize=fontsize_title)

#   Plot the colorbar

    if cbar_label != "" and abs(1-gpow) > 0.01:
        if not latex: gpow_label =   "^%4.2f"  %(gpow)
        else:         gpow_label = r'$^%4.2f$' %(gpow)
        cbar_label = cbar_label + gpow_label
        cbar_label_scale = 0.9
    else:
        cbar_label_scale = 1

    if plot_colorbar and plot_image:
        if   geometry == "annular":
            cbar = fig.colorbar(image,orientation='vertical',fraction=fraction_colorbar,pad=pad,aspect=aspect_cbar)
        elif geometry == "cartesian":
            cbar = plt.colorbar(image,orientation='vertical',fraction=fraction_colorbar,pad=pad,aspect=aspect_cbar)

        image.figure.axes[1].tick_params(axis="y",labelsize=fontsize_axis_tics)
        cbar.set_label(cbar_label,fontsize=cbar_label_scale*fontsize_axis_labels)
        if what in ("nu_kinematic","nu_dynamic"):
            import matplotlib.ticker as ticker
            # Force integer-only tick positions so each label is a unique power
            # of 10.  Without this, matplotlib places ticks at fractional log
            # values (e.g. 18.2, 18.5, 18.8) and the formatter rounds them all
            # to the same integer, producing duplicate "1e+XX" labels.
            cbar.ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
            def _visc_fmt(x, _):
                if not np.isfinite(x): return ""
                exp = int(round(x))
                return f"1e+{exp:02d}" if exp >= 0 else f"1e-{abs(exp):02d}"
            cbar.ax.yaxis.set_major_formatter(ticker.FuncFormatter(_visc_fmt))
            cbar.ax.set_title(cbar_label, fontsize=cbar_label_scale*fontsize_axis_labels, pad=4)
            cbar.set_label("")

#   Plot the axes labels and tick values

    if geometry == "cartesian":
        plt.xlim(x_min,x_max)
        plt.ylim(z_min, z_max)
        plt.tick_params(axis='x',labelsize=fontsize_axis_tics)
        plt.tick_params(axis='y',labelsize=fontsize_axis_tics)
       
#       Axes in lattice units

        if output_units == "lattice":
            if plane == "xz":
                if not latex:
                    plt.xlabel("x (lattice units)",fontsize=fontsize_axis_labels)
                    plt.ylabel("z (lattice units)",fontsize=fontsize_axis_labels)
                else:
                    plt.xlabel(r'$x$ (lattice units)',fontsize=fontsize_axis_labels)
                    plt.ylabel(r'$z$ (lattice units)',fontsize=fontsize_axis_labels)
            elif plane == "xy":
                if not latex:
                    plt.xlabel("x (lattice units)",fontsize=fontsize_axis_labels)
                    plt.ylabel("y (lattice units)",fontsize=fontsize_axis_labels)
                else:
                    plt.xlabel(r'$x$ (lattice units)',fontsize=fontsize_axis_labels)
                    plt.ylabel(r'$y$ (lattice units)',fontsize=fontsize_axis_labels)
            if plane == "yz":
                if not latex:
                    plt.xlabel("y (lattice units)",fontsize=fontsize_axis_labels)
                    plt.ylabel("z (lattice units)",fontsize=fontsize_axis_labels)
                else:
                    plt.xlabel(r'$y$ (lattice units)',fontsize=fontsize_axis_labels)
                    plt.ylabel(r'$z$ (lattice units)',fontsize=fontsize_axis_labels)

#       Axes in km x 1000 for large models such as the whole mantle

        elif L_p >= 1000000:
            if plane == "xz":
                if not latex:
                    plt.xlabel("x (km x 1000)",fontsize=fontsize_axis_labels)
                    plt.ylabel("z (km x 1000)",fontsize=fontsize_axis_labels)
                else:
                    plt.xlabel(r'$\displaystyle x$ (km \times 1000)',fontsize=fontsize_axis_labels)
                    plt.ylabel(r'$\displaystyle z$ (km \times 1000)',fontsize=fontsize_axis_labels)
            elif plane == "xy":
                if not latex:
                    plt.xlabel("x (km x 1000)",fontsize=fontsize_axis_labels)
                    plt.ylabel("y (km x 1000)",fontsize=fontsize_axis_labels)
                else:
                    plt.xlabel(r'$\displaystyle x$ (km \times 1000)',fontsize=fontsize_axis_labels)
                    plt.ylabel(r'$\displaystyle y$ (km \times 1000)',fontsize=fontsize_axis_labels)
            elif plane == "yz":
                if not latex:
                    plt.xlabel("y (km x 1000)",fontsize=fontsize_axis_labels)
                    plt.ylabel("z (km x 1000)",fontsize=fontsize_axis_labels)
                else:
                    plt.xlabel(r'$\displaystyle y$ (km \times 1000)',fontsize=fontsize_axis_labels)
                    plt.ylabel(r'$\displaystyle z$ (km \times 1000)',fontsize=fontsize_axis_labels)

#       Axes in km for small models such as the upper mantle

        else:
            if plane == "xz":
                if not latex:
                    plt.xlabel("x (km)",fontsize=fontsize_axis_labels)
                    plt.ylabel("z (km)",fontsize=fontsize_axis_labels)
                else:
                    plt.xlabel(r'$\displaystyle x$ (km)',fontsize=fontsize_axis_labels)
                    plt.ylabel(r'$\displaystyle z$ (km)',fontsize=fontsize_axis_labels)
            elif plane == "xy":
                if not latex:
                    plt.xlabel("x (km x 1000)",fontsize=fontsize_axis_labels)
                    plt.ylabel("y (km x 1000)",fontsize=fontsize_axis_labels)
                else:
                    plt.xlabel(r'$\displaystyle x$ (km)',fontsize=fontsize_axis_labels)
                    plt.ylabel(r'$\displaystyle y$ (km)',fontsize=fontsize_axis_labels)
            elif plane == "yz":
                if not latex:
                    plt.xlabel("y (km x 1000)",fontsize=fontsize_axis_labels)
                    plt.ylabel("z (km x 1000)",fontsize=fontsize_axis_labels)
                else:
                    plt.xlabel(r'$\displaystyle y$ (km)',fontsize=fontsize_axis_labels)
                    plt.ylabel(r'$\displaystyle z$ (km)',fontsize=fontsize_axis_labels)

#   Print the run description on the plot

    if plot_image == False: col_description = "k"
    else:                   col_description = description_color
    if description_location == "image":
        if description_box:
            plt.text(x_max-0.01*(x_max-x_min),z_min-0.14*(z_min-z_max),description,fontsize=fontsize_axis_labels,color=col_description,horizontalalignment='right',fontweight='bold',bbox=dict(boxstyle='square,pad=0.3',facecolor='white',edgecolor='black',linewidth=0.5))
        else:
            if col_description == "white" or col_description == "w": weight = 'black'
            else:                                                    weight = 'semibold'
            plt.text(x_max-0.01*(x_max-x_min),z_min-0.09*(z_min-z_max),description,fontsize=fontsize_axis_labels,color=col_description,horizontalalignment='right',fontweight=weight)
    elif description_location == "frame":
        if aspect <= 1.5: yshift = -0.03
        else:             yshift =  0.0
        if description_box:

            plt.figtext(0.99,0.06+yshift,description,fontsize=fontsize_axis_labels,color=col_description,horizontalalignment='right',fontweight='bold',bbox=dict(boxstyle='square,pad=0.3',facecolor='white',edgecolor='black',linewidth=0.5))
        else:
            if col_description == "white" or col_description == "w": weight = 'black'
            else:                                                    weight = 'semibold'
            plt.figtext(0.99,0.04+yshift,description,fontsize=fontsize_axis_labels,color=col_description,horizontalalignment='right',fontweight=weight)

#   Set the filename and save the figure

    if D == 2:
        what_file = what
    else:
        what_file = what + "_plane=" + plane 
        if   plane == "xz": what_file = what_file + "_y=" + str(coordinate)
        elif plane == "xy": what_file = what_file + "_z=" + str(coordinate)
        elif plane == "yz": what_file = what_file + "_x=" + str(coordinate)

    if t < 0: full_filename = plot_dir + what_file + ".png"
    else:     full_filename = plot_dir + what_file + "_%s.png" %(timetext(it_plots,timescale_Inc_plots,dt_p))

    try: os.remove(full_filename)
    except OSError: pass
    plt.savefig(full_filename, bbox_inches='tight')
    plt.close()
    if print_filename:
        print(", %s saved!" %(full_filename),flush=True)
    else:
        print()

#===============================================================================
# Utility routines
#===============================================================================

#-------------------------------------------------------------------------------
#   Print the range of a field (min and max)
#-------------------------------------------------------------------------------

def print_field_range(name,field):

    field_max,field_min = get_max_min_field(field)
    if rank == 0:
        print("%s E [%8.3g,%8.3g]" %(name,field_min,field_max),flush=True)

#-------------------------------------------------------------------------------
# Get max of a field and put onto core 0
#-------------------------------------------------------------------------------

def get_max_field(field):

    global Rx, Ry, Rz, mx, my, mz, xedge, yedge, zedge

    min_float = np.finfo(field.dtype).min

    if rank < nr:
        max_field = np.max(field[yedge:my[Ry]+yedge,zedge:mz[Rz]+zedge,xedge:mx[Rx]+xedge])*np.ones(1,dtype=field.dtype)
    else:
        max_field = min_float*np.ones(1,dtype=field.dtype)

    max_value = np.zeros(1,dtype=field.dtype)

    comm.Reduce(max_field,max_value,op=MPI.MAX,root=0)
    comm.Barrier()

    if rank == 0:
        return max_value[0]
    else:
        return 0

#-------------------------------------------------------------------------------
# Get max of a field including edges and put onto core 0
#-------------------------------------------------------------------------------

def get_max_field_all(field):

    min_float = np.finfo(field.dtype).min

    if rank < nr:
        max_field = np.max(field)*np.ones(1,dtype=field.dtype)
    else:
        max_field = min_float*np.ones(1,dtype=field.dtype)

    max_value = np.zeros(1,dtype=field.dtype)

    comm.Reduce(max_field,max_value,op=MPI.MAX,root=0)
    comm.Barrier()

    if rank == 0:
        return max_value[0]
    else:
        return 0

#-------------------------------------------------------------------------------
# Get min of a field and put onto core 0
#-------------------------------------------------------------------------------

def get_min_field(field):

    global Rx, Ry, Rz, mx, my, mz, xedge, yedge, zedge

    max_float = np.finfo(field.dtype).max

    if rank < nr:
        min_field = np.min(field[yedge:my[Ry]+yedge,zedge:mz[Rz]+zedge,xedge:mx[Rx]+xedge])*np.ones(1,dtype=field.dtype)
    else:
        min_field = max_float*np.ones(1,dtype=field.dtype)

    min_value = np.zeros(1,dtype=field.dtype)

    comm.Reduce(min_field,min_value,op=MPI.MAX,root=0)
    comm.Barrier()

    if rank == 0:
        return min_value[0]
    else:
        return 0

#-------------------------------------------------------------------------------
# Get min of a field including edges and put onto core 0
#-------------------------------------------------------------------------------

def get_min_field_all(field):

    max_float = np.finfo(field.dtype).max

    if rank < nr:
        min_field = np.min(field)*np.ones(1,dtype=field.dtype)
    else:
        min_field = max_float*np.ones(1,dtype=field.dtype)

    min_value = np.zeros(1,dtype=field.dtype)

    comm.Reduce(min_field,min_value,op=MPI.MAX,root=0)
    comm.Barrier()

    if rank == 0:
        return min_value[0]
    else:
        return 0

#-------------------------------------------------------------------------------
# Get max and min of a field and put onto core 0
#-------------------------------------------------------------------------------

def get_max_min_field(field):

    global Rx, Ry, Rz, mx, my, mz, xedge, yedge, zedge

    max_float = np.finfo(field.dtype).max
    min_float = np.finfo(field.dtype).min

    if rank < nr:
        min_field = np.min(field[yedge:my[Ry]+yedge,zedge:mz[Rz]+zedge,xedge:mx[Rx]+xedge])*np.ones(1,dtype=field.dtype)
        max_field = np.max(field[yedge:my[Ry]+yedge,zedge:mz[Rz]+zedge,xedge:mx[Rx]+xedge])*np.ones(1,dtype=field.dtype)
    else:
        min_field = max_float*np.ones(1,dtype=field.dtype)
        max_field = min_float*np.ones(1,dtype=field.dtype)

    min_value = np.zeros(1,dtype=field.dtype)
    max_value = np.zeros(1,dtype=field.dtype)

    comm.Reduce(max_field,max_value,op=MPI.MAX,root=0)
    comm.Reduce(min_field,min_value,op=MPI.MIN,root=0)
    comm.Barrier()

    if rank == 0:
        return max_value[0],min_value[0]
    else:
        return 0,0

#-------------------------------------------------------------------------------
# Get max and min of a field including edges and put onto core 0
#-------------------------------------------------------------------------------

def get_max_min_field_all(field):

    global Rx, Ry, Rz, mx, my, mz, xedge, yedge, zedge

    max_float = np.finfo(field.dtype).max
    min_float = np.finfo(field.dtype).min

    if rank < nr:
        min_field = np.min(field)*np.ones(1,dtype=field.dtype)
        max_field = np.max(field)*np.ones(1,dtype=field.dtype)
    else:
        min_field = max_float*np.ones(1,dtype=field.dtype)
        max_field = min_float*np.ones(1,dtype=field.dtype)

    min_value = np.zeros(1,dtype=field.dtype)
    max_value = np.zeros(1,dtype=field.dtype)

    comm.Reduce(max_field,max_value,op=MPI.MAX,root=0)
    comm.Reduce(min_field,min_value,op=MPI.MIN,root=0)
    comm.Barrier()

    if rank == 0:
        return max_value[0],min_value[0]
    else:
        return 0,0

#-------------------------------------------------------------------------------
# Get the mean of a field and put onto core 0
#-------------------------------------------------------------------------------

def get_mean_field(field):

    global Rx, Ry, Rz, mx, my, mz, xedge, yedge, zedge

    if rank < nr:
        sum_field  = np.sum(field[yedge:my[Ry]+yedge,zedge:mz[Rz]+zedge,xedge:mx[Rx]+xedge])*np.ones(1,dtype=field.dtype)
    else:
        sum_field = np.zeros(1,dtype=field.dtype)

    mean_field = np.zeros(1,dtype=field.dtype)

    comm.Reduce(sum_field,mean_field,op=MPI.SUM,root=0)
    comm.Barrier()

    if rank == 0:
        return mean_field[0]/(nx*ny*nz)
    else:
        return 0

#-------------------------------------------------------------------------------
# Get max value on used cores to rank 0
#-------------------------------------------------------------------------------

def get_max_value(value_scalar):

    value    = np.zeros(1,dtype=np.float32)
    max_float = np.finfo(value.dtype).max

    if rank < nr:
        value[0] = value_scalar
    else:
        value[0] = - max_float

    max_value = np.zeros(1,dtype=np.float32)

    comm.Reduce(value,max_value,op=MPI.MAX,root=0)
    comm.Barrier()

    if rank == 0:
        return max_value[0]
    else:
        return 0

#-------------------------------------------------------------------------------
# Get mean value on used cores to rank 0
#-------------------------------------------------------------------------------

def get_mean_value(value):

    if rank == 0:
        values    = np.zeros(nr,np.float32)
        values[0] = value
    elif rank < nr: comm.send(value,dest=0)

    if rank == 0:
        for r in range(1,nr): values[r] = comm.recv(source=r)

        values_sum = np.sum(values)
        mean       = values_sum/nr

        return mean
    else:
        return 0

#-------------------------------------------------------------------------------
# Print performance and timing
#-------------------------------------------------------------------------------

def print_timing():

    global totalTimeEdges, totalTimeEdges1, totalTimeEdges2, totalTimeStream, totalTimeCalcs, totalTimeCalcs1, totalTimeCalcs2, totalTimeMach, totalTimeOutput, time_h, time_s, time_m, Mflups_LBM, Mflups_elapsed

    elapsed_time        = time.time() - time0

    elapsed_time        = get_max_value (elapsed_time)

    totalTimeEdges   = get_max_value(totalTimeEdges)
    totalTimeEdges1  = get_max_value(totalTimeEdges1)
    totalTimeEdges2  = get_max_value(totalTimeEdges2)
    totalTimeCalcs   = get_max_value(totalTimeCalcs)
    totalTimeCalcs1  = get_max_value(totalTimeCalcs1)
    totalTimeCalcs2  = get_max_value(totalTimeCalcs2)
    totalTimeStream  = get_max_value(totalTimeStream)
    totalTimeMach    = get_max_value(totalTimeMach)
    totalTimeOutput  = get_max_value(totalTimeOutput)

    if rank == 0:

        elapsed_1000        = elapsed_time /((t+1)/1000)

        totalTimeLBM     = totalTimeEdges  + totalTimeCalcs + totalTimeMach
        totalTime        = totalTimeLBM    + totalTimeOutput

#       Correct for drift because adding many small values to get timing breakdown

        correction           = (elapsed_time+zero)/(totalTime+zero)

        totalTimeEdges   *= correction
        totalTimeEdges1  *= correction
        totalTimeEdges2  *= correction
        totalTimeCalcs   *= correction
        totalTimeCalcs1  *= correction
        totalTimeCalcs2  *= correction
        totalTimeStream  *= correction
        totalTimeMach    *= correction
        totalTimeOutput  *= correction
        totalTimeLBM     *= correction
        totalTime        *= correction

        totalTime_1000       = totalTime /((t+1)/1000)

        flups_per_core = float(NX)*(t_computed+1)
        flops_per_core = float(Flops)*flups_per_core

        Mflups_LBM     = flups_per_core/(1.e6*totalTimeLBM)
        Mflups_elapsed = flups_per_core/(1.e6*elapsed_time)
        Gflops_LBM     = flops_per_core/(1.e9*totalTimeLBM)
        Gflops_elapsed = flops_per_core/(1.e9*elapsed_time)

        time_s = elapsed_time
        time_m = int(time_s/60)
        time_h = int(time_m/60)
        time_m = time_m - time_h*60
        time_s = time_s - time_m*60 - time_h*60*60

        now = datetime.now()
        date_string =  now.strftime("%B %d, %Y, %H:%M:%S")

        if it_timing == nt:
            if it_timing != it_output:
                print("==================================================================",flush=True)
            print("FULL RUN TIMING/PERFORMANCE, run/project = %s/%s at %s" %(RUN,PROJECT,date_string),flush=True)
            print("==================================================================",flush=True)
            print("Cores specified              = ",size,flush=True)
            print("Cores used                   = ",nr,flush=True)
            print("------------------------------------------------------------------",flush=True)
            print("domain decomposition         = ",domain_decomposition,flush=True)
            print("------------------------------------------------------------------",flush=True)
            print("code (python | c)            = ",code,flush=True)
            print("overlap                      = ",overlap,flush=True)
            print("float size calc              = ",float_bits,flush=True)
            print("float size edges             = ",edge_bits,flush=True)
            print("save size                    = ",data_bits,flush=True)
            print("------------------------------------------------------------------",flush=True)
            print("Model  size                  = (",nx," x ",ny," x ",nz,")",flush=True)
            print("Domain size                  = (",nxi," x ",nyi," x ",nzi,")",flush=True)
            print("Domains                      = (",nrx," x ",nry," x ",nrz,")",flush=True)
            print("------------------------------------------------------------------",flush=True)
            print("Domain size in lattice sites = ",NX,flush=True)
            print("Flops/timestep/lattice site  = ",Flops,flush=True)
        else:
            print("==================================================================",flush=True)
            print("TIMING & PERFORMANCE of run/project = %s/%s at %s" %(RUN,PROJECT,date_string),flush=True)

        print("==================================================================",flush=True)
        print("Timing breakdown",flush=True)
        print("==================================================================",flush=True)

        if overlap and code == "c":
            print("Edges start:    %9.1f secs"%(totalTimeEdges1),   " = %4.1f"%(100*totalTimeEdges1/totalTimeEdges),"% of edges",flush=True)
            print("Edges end:      %9.1f secs"%(totalTimeEdges2),   " = %4.1f"%(100*totalTimeEdges2/totalTimeEdges),"% of edges",flush=True)
            print("------------------------------------------------------------------",flush=True)
            print("Calcs interior: %9.1f secs"%(totalTimeCalcs1),  " = %4.1f"%(100*totalTimeCalcs1/totalTimeCalcs),"% of calcs",flush=True)
            print("Calcs edges:    %9.1f secs"%(totalTimeCalcs2),  " = %4.1f"%(100*totalTimeCalcs2/totalTimeCalcs),"% of calcs",flush=True)
            print("==================================================================",flush=True)
        print("Edges:          %9.1f secs"%(totalTimeEdges),   " = %4.1f"%(100*totalTimeEdges/totalTimeLBM),"% of LBM",flush=True)
        print("Calcs:          %9.1f secs"%(totalTimeCalcs),  " = %4.1f"%(100*totalTimeCalcs/totalTimeLBM),"% of LBM",flush=True)
        print("Mach:           %9.1f secs"%(totalTimeMach),   " = %4.1f"%(100*totalTimeMach/totalTimeLBM),"% of LBM",flush=True)
        if output_plots or output_data:
            print("------------------------------------------------------------------",flush=True)
            print("Edges:          %9.1f secs"%(totalTimeEdges),  " = %4.1f"%(100*totalTimeEdges/totalTime),"% of total",flush=True)
            print("Calcs:          %9.1f secs"%(totalTimeCalcs),  " = %4.1f"%(100*totalTimeCalcs/totalTime),"% of total",flush=True)
            print("Mach:           %9.1f secs"%(totalTimeMach),   " = %4.1f"%(100*totalTimeMach/totalTime),"% of total",flush=True)
            print("Output:         %9.1f secs"%(totalTimeOutput),  " = %4.1f"%(100*totalTimeOutput/totalTime),"% of total",flush=True)
        else:
            print("------------------------------------------------------------------",flush=True)
            print("Mach:           %9.1f secs"%(totalTimeMach),   " = %4.1f"%(100*totalTimeMach/totalTime),"% of total",flush=True)
        print("==================================================================",flush=True)
        print("Performance: ncores = ",nr,flush=True)
        print("==================================================================",flush=True)
        print("Gflops/core (excluding output)   = %6.2f" %(Gflops_LBM),flush=True)
        print("Gflops/core                      = %6.2f" %(Gflops_elapsed),flush=True)
        print("Gflops                           = %6.2g" %(Gflops_elapsed*nr),flush=True)
        print("------------------------------------------------------------------",flush=True)
        print("Mflups/core (excluding output)   = %6.2f" %(Mflups_LBM),flush=True)
        print("Mflups/core                      = %6.2f" %(Mflups_elapsed),flush=True)
        print("Mflups                           = %6.2g" %(Mflups_elapsed*nr),flush=True)
        if it_timing < nt:
            print("------------------------------------------------------------------",flush=True)
            print("Time                             = ",it_timing,flush=True)
            if output_units == "physical":
                time_text = timesteps_to_physical_time_text(t,dt_p,output_units)
                print("Time = %d (%d computed) = %s" %(it_mach,t_computed,time_text),flush=True)
            else:
                print("Time = %d (%d computed)" %(it_mach,t_computed),flush=True)
            print("------------------------------------------------------------------",flush=True)
            print("Secs/Ktimestep                   =  %5.3f secs" %(elapsed_1000),flush=True)
            print("Run time: %d hrs, %d mins, %4.1f secs" %(time_h,time_m,time_s),flush=True)
        elif it_timing == nt:
            for iout in range(2):
                if   iout == 0: output_print = sys.stderr
                elif iout == 1: output_print = sys.stdout
                print("==================================================================",file=output_print,flush=True)
                print("Run complete: run/project = %s/%s at %s" %(RUN,PROJECT,date_string),file=output_print,flush=True)
                print("------------------------------------------------------------------",flush=True)
                print("Accuracy health check: mach(ave,max) = (%5.3f,%5.3f)" %(mach_max_ave,mach_max_run),end="",flush=True)
                if mach_max_ave > 0.05: print(", Warning for high Ma: Ma_ave > 0.05",end="",flush=True)
                if mach_max_run > 0.05: print(", Warning for high Ma: Ma_max > 0.05",end="",flush=True)
                print("",flush=True)
                print("==================================================================",file=output_print,flush=True)
                print("Run time: %d hrs, %d mins, %4.1f secs" %(time_h,time_m,time_s),file=output_print,flush=True)
        print("==================================================================",flush=True)

#===============================================================================
# Initialize distributions
#===============================================================================

#-------------------------------------------------------------------------------
# Initialize distributions fast on cores by assuming u = 0
#-------------------------------------------------------------------------------

def initialize_distributions_fast():
    if RKTLBM :
       u2  = np.einsum('yzdx,yzdx->yzx',u,u)
       rho = rho_r/rho1 + rho_b/rho2
       for a in range(na):
              cu = np.einsum('d,yzdx->yzx',c[a],u)
              if shifted_distribution:
                  if RKTLBM:
                      # Shifted RKTLBM: reference must be the spatially UNIFORM constant rho1/rho2
                      # (set in rho_r_init = np.full(..., rho1) below) so that
                      # rho_r_init[source] = rho_r_init[dest] during streaming and the MAC error
                      # term (rho_r_init[dest] − Σ_a rho_r_init[src_a]·w[a]) = 0 exactly.
                      # Initialise to the full shifted equilibrium: f̃_r = feq_r − rho1·C_r[a].
                      f_r [:,:,a,:] = rho_r  * (C_r[a] + w[a] * (    3*cu + 4.5*cu**2 - 1.5*u2)) - rho1 * C_r[a]
                      f_b [:,:,a,:] = rho_b  * (C_b[a] + w[a] * (    3*cu + 4.5*cu**2 - 1.5*u2)) - rho2 * C_b[a]
                  else:
                      # Single-fluid shifted: reference ≈ 1 (rho≈1), f̃ ≈ 0 at u=0.
                      f_r [:,:,a,:] = rho_r  * w[a] * (    3*cu + 4.5*cu**2 - 1.5*u2)
                      f_b [:,:,a,:] = rho_b  * w[a] * (    3*cu + 4.5*cu**2 - 1.5*u2)
              else:
                  f_r [:,:,a,:] = rho_r  * (C_r[a] + w[a] * (    3*cu + 4.5*cu**2 - 1.5*u2))
                  f_b [:,:,a,:] = rho_b  * (C_b[a] + w[a] * (    3*cu + 4.5*cu**2 - 1.5*u2))
              g   [:,:,a,:] = eta    * rho     * w[a] * (1 + 3*cu + 4.5*cu**2 - 1.5*u2)
              if rheology == "archean2":
                 h[:,:,a,:] = L_ed             * w[a] * (1 + 3*cu + 4.5*cu**2 - 1.5*u2)
    else:
       for a in range(na):
           if shifted_distribution: f[:,:,a,:] = 0
           else:                    f[:,:,a,:] = w[a]
           g[:,:,a,:] = eta * w[a]
           if rheology == "archean2":
              h[:,:,a,:] = L_ed * w[a]
            
#-------------------------------------------------------------------------------
# Initialize distributions on cores for any initial velocity u
#-------------------------------------------------------------------------------

def initialize_distributions():

    u2 = np.einsum('yzdx,yzdx->yzx',u,u)
    for a in range(na):
        if D == 2:
            cu = np.einsum('d,yzdx->yzx',c[a],u)
        else:
            cu = np.einsum('d,yzdx->yzx',c[a],u)

        if RKTLBM :
              rho = rho_r/rho1+rho_b/rho2
              f_r[:,:,a,:]  = rho_r * (C_r[a] + w[a] * (3*cu + 4.5*cu**2 - 1.5*u2))
              f_b[:,:,a,:]  = rho_b * (C_b[a] + w[a] * (3*cu + 4.5*cu**2 - 1.5*u2))
              g[:,:,a,:]    = eta   * rho     * w[a] * (1 + 3*cu + 4.5*cu**2 - 1.5*u2)
              if rheology == "archean2":
                 h[:,:,a,:] = L_ed  * w[a] * (1 + 3*cu + 4.5*cu**2 - 1.5*u2)

        else :  
           if shifted_distribution:
              f[:,:,a,:] = rho * w[a] * (3*cu + 4.5*cu**2 - 1.5*u2) + w[a] * (rho-1)
              g[:,:,a,:] = eta * (f[:,:,a,:] + rho*w[a])
           else:
              f[:,:,a,:] = rho * w[a] * (1 + 3*cu + 4.5*cu**2 - 1.5*u2)
              g[:,:,a,:] = eta * f[:,:,a,:]
              if rheology == "archean2":
                 h[:,:,a,:] = L_ed * w[a] * (1 + 3*cu + 4.5*cu**2 - 1.5*u2)

#===============================================================================
# Domain decomposition dependent routines
#===============================================================================

import optimal_domain_decomposition
from optimal_domain_decomposition import my_domain_decomposition

#-------------------------------------------------------------------------------
# Initialize fast get_field which uses persistent send/recv requests
#-------------------------------------------------------------------------------

def initialize_get_field_fast():

    global reqSend_get_field, reqRecv_get_field, get_field_send, get_field_recv
    global nr, nxi, nyi, nzi

    if rank == 0:
        get_field_recv = np.zeros((nr-1,nyi,nzi,nxi),dtype=np_data) # Recv field on rank 0 for rank > 0
    elif rank < nr:
        get_field_send = np.zeros(     (nyi,nzi,nxi),dtype=np_data) # Send field on ranks > 0 to core 0

#   Initialize persistent send and recv requests

    if rank > 0 and rank < nr:
        reqSend_get_field       = comm.Send_init([get_field_send   ,MPI_DATA],dest  =0,tag=rank)
    elif rank == 0:
        reqRecv_get_field = ["NULL"]
        for r in range(1,nr):
            reqRecv_get_field_r = comm.Recv_init([get_field_recv[r-1],MPI_DATA],source=r,tag=r)
            reqRecv_get_field.append(reqRecv_get_field_r)

#-------------------------------------------------------------------------------
# Initialize fast put_field using persistent send/recv requests
#-------------------------------------------------------------------------------

def initialize_put_field_fast():

    global reqSend_put_field, reqRecv_put_field, put_field_send, put_field_recv
    global nr, nxi, nyi, nzi

    if rank == 0:
        put_field_send = np.zeros((nr-1,nyi,nzi,nxi),dtype=np_data) # Send field on rank 0 to ranks > 0
    elif rank < nr:
        put_field_recv = np.zeros(     (nyi,nzi,nxi),dtype=np_data) # Recv field on ranks > 0 from core 0

#   Initialize persistent send and recv requests

    if rank > 0 and rank < nr:
        reqRecv_put_field       = comm.Recv_init([put_field_recv     ,MPI_DATA],source  =0,tag=rank)
    elif rank == 0:
        reqSend_put_field = ["NULL"]
        for r in range(1,nr):
            reqSend_put_field_r = comm.Send_init([put_field_send[r-1],MPI_DATA],dest    =r,tag=r)
            reqSend_put_field.append(reqSend_put_field_r)

#-------------------------------------------------------------------------------
# Fast get_field to core 0 using persistent send/recv requests
#-------------------------------------------------------------------------------

def get_field(field):

    global reqSend_get_field, reqRecv_get_field, get_field_send, get_field_recv
    global xedge, yedge, zedge, mx, my, mz, nxi, nyi, nzi

#   Before sending, convert float_bits precision to fit within the precision
#   of data_bits precision without going to infinity which will mess up plots

    data_smallest = np.finfo(field.dtype).min
    data_largest  = np.finfo(field.dtype).max
    field         = np.where(field<data_smallest,data_smallest,field)
    field         = np.where(field>data_largest ,data_largest,field )

#   Start all sends and receives

    if rank > 0 and rank < nr:
        get_field_send[:,:,:] = field[yedge:yedge+nyi,zedge:zedge+nzi,xedge:xedge+nxi]
        MPI.Prequest.Startall([reqSend_get_field])
    if rank == 0:
        for r in range(1,nr):
            MPI.Prequest.Startall([reqRecv_get_field[r]])

#   Wait for all sends and receives to finish

    if rank > 0 and rank < nr:
        MPI.Request.Waitall     ([reqSend_get_field])
    elif rank == 0:
        for r in range(1,nr):
            MPI.Request.Waitall([reqRecv_get_field[r]])

#   Put each domain into the correct positon in the FIELD array

    if rank == 0:
        FIELD      = np.zeros((ny,nz,nx),dtype=np_data)
        FIELD[0:nyi,0:nzi,0:nxi] = field[yedge:yedge+nyi,zedge:zedge+nzi,xedge:xedge+nxi]
        for rx in range(nrx):
            for ry in range(nry):
                for rz in range(nrz):
                    r = rx + rz*nrx + ry*nrx*nrz
                    x0 = rx * nxi
                    x1 = x0 + mx[rx]
                    y0 = ry * nyi
                    y1 = y0 + my[ry]
                    z0 = rz * nzi
                    z1 = z0 + mz[rz]
                    if r > 0:
                        FIELD[y0:y1,z0:z1,x0:x1] = get_field_recv[r-1,0:my[ry],0:mz[rz],0:mx[rx]]
    else:
        FIELD = 0


    return FIELD

#-------------------------------------------------------------------------------
# Fast put_field from core 0 using persistent send/recv requests
#
# Note: The default value is put into the edges and does not affect calculations
# which are only done in the physical space for code = c. However, for
# code = python, and considering python's numpy operates on the entire
# array including edges, a value that avoids python overflows should be
# specified as the default value.
#-------------------------------------------------------------------------------

def put_field(FIELD,default=0):

    global reqSend_put_field, reqRecv_put_field, put_field_send, put_field_recv
    global xedge, yedge, zedge, mx, my, mz, nxi, nyi, nzi
    global Rx, Ry, Rz

#   Initialize field on domains

    if rank < nr: field = default * np.ones((ny_i,nz_i,nx_i),dtype=np_float)

#   Start all sends and receives

    if rank == 0:
        for rx in range(nrx):
            for ry in range(nry):
                for rz in range(nrz):
                    r = rx + rz*nrx + ry*nrx*nrz
                    x0 = rx * nxi
                    x1 = x0 + mx[rx]
                    y0 = ry * nyi
                    y1 = y0 + my[ry]
                    z0 = rz * nzi
                    z1 = z0 + mz[rz]
                    if r > 0:
                        put_field_send[r-1,:my[ry],:mz[rz],:mx[rx]] = FIELD[y0:y1,z0:z1,x0:x1]
                        MPI.Prequest.Startall([reqSend_put_field[r]])
    elif rank < nr:
        MPI.Prequest.Startall([reqRecv_put_field])

#   Wait for all sends and receives to finish

    if rank == 0:
        for r in range(1,nr):
            MPI.Request.Waitall([reqSend_put_field[r]])
    elif rank < nr:
        MPI.Request.Waitall([reqRecv_put_field])

#   Put each domain into the correct positon in the FIELD array

    if rank == 0:
        field[yedge:yedge+nyi,zedge:zedge+nzi,xedge:xedge+nxi] = FIELD[0:nyi,0:nzi,0:nxi]
    elif rank < nr:
        field         [yedge:yedge+my[Ry],zedge:zedge+mz[Rz],xedge:xedge+mx[Rx]] = \
        put_field_recv[0:my[Ry]          ,0:mz[Rz]          ,0:mx[Rx]]
    else:
        field = 0

    return field

#-------------------------------------------------------------------------------
# get_field including the edges onto core 0 (used for testing/debugging)
#-------------------------------------------------------------------------------

def get_field_all(field):

    global R, Rx, Ry, Rz

    if rank != 0 and rank < nr:
        comm.send(field,dest=0)

    if rank == 0:
        FIELD = np.zeros((nry*ny_i,nrz*nz_i,nrx*nx_i),dtype=field.dtype)
        for rx in range(nrx):
            for ry in range(nry):
                for rz in range(nrz):
                    r = rx + rz*nrx + ry*nrx*nrz
                    x0 = rx * nx_i
                    x1 = x0 + nx_i
                    y0 = ry * ny_i
                    y1 = y0 + ny_i
                    z0 = rz * nz_i
                    z1 = z0 + nz_i
                    if r == 0:
                        FIELD[y0:y1,z0:z1,x0:x1] = field
                    else:
                        tmp = comm.recv(source=r)
                        FIELD[y0:y1,z0:z1,x0:x1] = tmp
    else:
        FIELD = False

    return FIELD

#-------------------------------------------------------------------------------
# Initialize sinusoidal temperature variations for benchmarking
#-------------------------------------------------------------------------------
 
def initialize_sinusoidal_temperature(T):

    global Rx, Ry, Rz

    x     = np.linspace(-xedge+Rx*nxi,(Rx+1)*nxi,nx_i,dtype=np_float)
    y     = np.linspace(-yedge+Ry*nyi,(Ry+1)*nyi,ny_i,dtype=np_float)
    z     = np.linspace(-zedge+Rz*nzi,(Rz+1)*nzi,nz_i,dtype=np_float)
    y,z,x = np.meshgrid(y,z,x,indexing="ij")

    L_x       = max(1,nx-1)
    L_y       = max(1,ny-1)
    a_x       = 2*np.pi*x_periods/L_x
    a_y       = 2*np.pi*y_periods/L_y
    a_z       = np.pi/L_l

    T = np_float(dT_p * ((z-L_l/2)/L_l + amp_sinusoid*np.sin(a_z*z)*np.cos(a_x*x)*np.cos(a_y*y)))

    return T

#-------------------------------------------------------------------------------
# Put a function of z on core 0 to a field which is a function of x, y and z
# on all cores
#
# Note: The default value is put into the edges and does not affect calculations
# which are only done in the physical space for code = c. However, for
# code = python, and considering python's numpy operates on the entire
# array including edges, a value that avoids python overflows should be
# specified as the default value.
#-------------------------------------------------------------------------------
 
def put_field_z_2_xyz(FIELD_z,default=0):

    global Rz

    field = default*np.ones((ny_i,nz_i,nx_i),dtype=np_float)
    for z in range(nz):
        rz = z//nzi
        if rz == Rz:
            z0  = rz * nzi
            z1  = z0 + mz[rz]

#          This loop needed, cannot broadcast FIELD_z[:] to field[:,:,:]

            for zi in range(z0,z1):
                iz = zedge+zi%nzi
                field[:,iz,:] = FIELD_z[zi]
                if xedge > 0:
                    field[:,iz,0       ] = FIELD_z[zi]
                    field[:,iz,mx[Rx]+1] = FIELD_z[zi]
                if yedge > 0:
                    field[0       ,iz,:] = FIELD_z[zi]
                    field[my[Ry]+1,iz,:] = FIELD_z[zi]

    return field

#-------------------------------------------------------------------------------
# Put a function of z on core 0 to a field which is a function of z only
# on all cores.

# Notes: (1) This is used for code = c which does not send a full
# array of xyz to the c code if the variable is only a function of z
#
#        (2) The default value is put into the edges and does not affect calculations
# which are only done in the physical space for code = c. However, for
# code = python, and considering python's numpy operates on the entire
# array including edges, a value that avoids python overflows should be
# specified as the default value.
#-------------------------------------------------------------------------------

def put_field_z_2_z(FIELD_z,default=0):

    global Rz

    field = default*np.ones((nz_i),dtype=np_float)
    for rx in range(nrx):
        for ry in range(nry):
            for rz in range(nrz):
                r = rx + rz*nrx + ry*nrx*nrz
                if rank == r and rz == Rz:
                    z0 = rz * nzi
                    z1 = z0 + mz[rz]
                    # Fill INNER cells from the global depth profile.
                    field[zedge:zedge+mz[rz]] = FIELD_z[z0:z1]
                    # Fill HALO cells from the neighboring slice of the
                    # global profile so the per-z values are continuous
                    # across rank boundaries (otherwise a multi-rank z
                    # decomposition leaves the halo at `default`, which
                    # produces a jump at the rank boundary in nu_f_b /
                    # Delta_eta_array / act_E / act_V — visible in the
                    # nu_kinematic npy as a discontinuity at z = nzi).
                    nz_global = FIELD_z.shape[0]
                    for k in range(zedge):
                        zsrc = z0 - 1 - k
                        if   zsrc >= 0:        field[zedge - 1 - k]              = FIELD_z[zsrc]
                        else:                  field[zedge - 1 - k]              = FIELD_z[0]
                        zsrc = z1 + k
                        if   zsrc < nz_global: field[zedge + mz[rz] + k]         = FIELD_z[zsrc]
                        else:                  field[zedge + mz[rz] + k]         = FIELD_z[nz_global - 1]
    return field

#-------------------------------------------------------------------------------
# Initialize BC
#-------------------------------------------------------------------------------

def initialize_BC():

    global xl        , xr        , yb        , yf        , zu        , zd
    global boundary_l, boundary_r, boundary_b, boundary_f, boundary_u, boundary_d
    global boundary_bb, boundary_ff, boundary_uu, boundary_dd
    global nbL, nbR, nbB, nbF, nbU, nbD
    global bL , bR , bB , bF , bU , bD
    global biL, biR, biB, biF, biU, biD
    global eta_u, eta_d
    global R, Rx, Ry, Rz

#   Set left and right BC variables

    xl         = 0
    xr         = nx_i - 1
    boundary_l = False
    boundary_r = False
    bL         = np.zeros(na,dtype=np_int)
    bR         = np.zeros(na,dtype=np_int)
    biL        = np.zeros(na,dtype=np_int)
    biR        = np.zeros(na,dtype=np_int)
    nbL        = naL
    nbR        = naR
    bL[:]      = aL[:]
    bR[:]      = aR[:]
    if BC_sides == "fixed" or BC_sides == "free_slip":
        if Rx == 0:
            xl         = xedge+0     %nxi
            boundary_l = True
        if Rx == nrx-1:
            xr          = xedge+(nx-1)%nxi
            boundary_r = True
        if   BC_sides == "fixed":
            for a in range(naL):
                biR[a] = ai[aR[a]]
                biL[a] = ai[aL[a]]
            for a in range(1,na):
                for b in range(1,na):
                    if c[a,X] == 0 and c[b,X] == 0:
                        if c[a,Z] == -c[b,Z] and (D == 2 or (D == 3 and c[a,Y] == -c[b,Y])):
                            if c[a,Z] >= 0 and (D == 2 or (D == 3 and c[a,Y] >= 0)):
                                bR [nbR] = a
                                biR[nbR] = b
                                nbR     += 1
                                bL [nbL] = a
                                biL[nbL] = b
                                nbL     += 1
        elif BC_sides == "free_slip":
            for a in range(naR):
                for b in range(naL):
                    if c[aR[a],X] == -c[aL[b],X]:
                        if (D == 2 and c[aR[a],Z] == c[aL[b],Z]) or (D == 3 and c[aR[a],Z] == c[aL[b],Z] and c[aR[a],Y] == c[aL[b],Y]):
                            biR[a] = bL[b]
                            biL[b] = bR[a]

#   Set backward and forward BC variables

    yb          = 0
    yf          = ny_i - 1
    boundary_b  = False
    boundary_f  = False
    boundary_bb = False
    boundary_ff = False
    bB         = np.zeros(na,dtype=np_int)
    bF         = np.zeros(na,dtype=np_int)
    biB        = np.zeros(na,dtype=np_int)
    biF        = np.zeros(na,dtype=np_int)
    nbB        = naB
    nbF        = naF
    bB[:]      = aB[:]
    bF[:]      = aF[:]
    if (BC_sides == "fixed" or BC_sides == "free_slip") and D == 3:
        if Ry == 0:
            yb         = yedge+0     %nyi
            boundary_b = True
        if Ry == nry-1:
            yf         = yedge+(ny-1)%nyi
            boundary_f = True

#       Backward and forward fixed or free_slip BC

        if   BC_sides == "fixed":
            for a in range(naF):
                biF[a] = ai[aF[a]]
                biB[a] = ai[aB[a]]
            for a in range(1,na):
                for b in range(1,na):
                    if c[a,Y] == 0 and c[b,Y] == 0:
                        if c[a,X] == -c[b,X] and c[a,Z] == -c[b,Z]:
                            if c[a,X] >= 0 and c[a,Z] >= 0:
                                bF [nbF] = a
                                biF[nbF] = b
                                nbF     += 1
                                bB [nbB] = a
                                biB[nbB] = b
                                nbB     += 1
        elif BC_sides == "free_slip":
            for a in range(naF):
                for b in range(naB):
                    if c[aF[a],Y] == -c[aB[b],Y]:
                        if c[aF[a],X] == c[aB[b],X] and c[aF[a],Z] == c[aB[b],Z]:
                            aiF[a] = aB[b]
                            aiB[b] = aF[a]

#   Set upper and lower BC variables

    zu          = 0
    zd          = nz_i - 1
    boundary_u  = False
    boundary_d  = False
    boundary_uu = False
    boundary_dd = False
    bU         = np.zeros(na,dtype=np_int)
    bD         = np.zeros(na,dtype=np_int)
    biU        = np.zeros(na,dtype=np_int)
    biD        = np.zeros(na,dtype=np_int)
    nbU        = naU
    nbD        = naD
    bU[:]      = aU[:]
    bD[:]      = aD[:]
    eta_u       = np.zeros((ny_i,nx_i),dtype=np_float) # Energy density upper boundary
    eta_d       = np.zeros((ny_i,nx_i),dtype=np_float) # Energy density lower boundary
    if Rz == 0:
        zu         = zedge+0     %nzi
        boundary_u = True
    if Rz == nrz-1:
        zd         = zedge+(nz-1)%nzi
        boundary_d = True

#   Upper fixed or free_slip BC

    if   BC_upper == "fixed":
        for a in range(naD):
            biD[a] = ai[bD[a]]
        for a in range(1,na):
            for b in range(1,na):
                if c[a,Z] == 0 and c[b,Z] == 0:
                    if c[a,X] == -c[b,X] and (D == 2 or (D == 3 and c[a,Y] == -c[b,Y])):
                        if c[a,X] >= 0 and (D == 2 or (D == 3 and c[a,Y] >= 0)):
                            bD [nbD] = a
                            biD[nbD] = b
                            nbD     += 1
    elif BC_upper == "free_slip":
        for a in range(naD):
            for b in range(naU):
                if c[aD[a],Z] == -c[aU[b],Z]:
                    if (D == 2 and c[aD[a],X] == c[aU[b],X]) or (D == 3 and c[aD[a],X] == c[aU[b],X] and c[aD[a],Y] == c[aU[b],Y]):
                        biD[a] = aU[b]

#   Lower fixed or free_slip BC

    if   BC_lower == "fixed":
        for a in range(naU):
            biU[a] = ai[bU[a]]
        for a in range(1,na):
            for b in range(1,na):
                if c[a,Z] == 0 and c[b,Z] == 0:
                    if c[a,X] == -c[b,X] and (D == 2 or (D == 3 and c[a,Y] == -c[b,Y])):
                        if c[a,X] >= 0 and (D == 2 or (D == 3 and c[a,Y] >= 0)):
                            bU [nbU] = a
                            biU[nbU] = b
                            nbU     += 1
    elif BC_lower == "free_slip":
        for a in range(naU):
            for b in range(naD):
                if c[aD[a],Z] == -c[aU[b],Z]:
                    if (D == 2 and c[aD[a],X] == c[aU[b],X]) or (D == 3 and c[aD[a],X] == c[aU[b],X] and c[aD[a],Y] == c[aU[b],Y]):
                        biU[a] = aD[b]

#    elif BC_lower == "free_slip":
#        for a in range(naU):
#            for b in range(naD):
#                if c[aU[a],Z] == -c[aD[b],Z]:
#                    if (D == 2 and c[aU[a],X] == c[aD[b],X]) or (D == 3 and c[aU[a],X] == c[aD[b],X] and c[aU[a],Y] == c[aD[b],Y]):
#                        biU[a] = aD[b]



#-------------------------------------------------------------------------------
# Apply fixed or free_slip BC as required (used for code = python runs)
#-------------------------------------------------------------------------------

def apply_BC():

    global xl        , xr        , yb        , yf        , zu        , zd
    global boundary_l, boundary_r, boundary_b, boundary_f, boundary_u, boundary_d
    global nbL, nbR, nbB, nbF, nbU, nbD
    global bL , bR , bB , bF , bU , bD
    global biL, biR, biB, biF, biU, biD

#   Left and right BC

    if boundary_l:

       if RKTLBM:
          for a in range(naR):
              f_r[:,:,bR[a],xl] = f_r[:,:,biR[a],xl]   
              f_b[:,:,bR[a],xl] = f_b[:,:,biR[a],xl]
              g[:,:,bR[a],xl]   = g  [:,:,biR[a],xl]
              if rheology == "archean2":
                 h[:,:,bR[a],xl] = h[:,:,biR[a],xl]
            
          for a in range(naR,nbR):
              tmp              = 0.5*(f_r[:,:,bR[a],xl] + f_r[:,:,biR[a],xl])
              f_r[:,:,bR [a],xl] = tmp
              f_r[:,:,biR[a],xl] = tmp
              tmp              = 0.5*(f_b[:,:,bR[a],xl] + f_b[:,:,biR[a],xl])
              f_b[:,:,bR [a],xl] = tmp
              f_b[:,:,biR[a],xl] = tmp
              tmp              = 0.5*(g[:,:,bR[a],xl] + g[:,:,biR[a],xl])
              g[:,:,bR [a],xl] = tmp
              g[:,:,biR[a],xl] = tmp
              if rheology == "archean2":
                 tmp              = 0.5*(h[:,:,bR[a],xl] + h[:,:,biR[a],xl])
                 h[:,:,bR [a],xl] = tmp
                 h[:,:,biR[a],xl] = tmp


       else:
          for a in range(naR):
              f[:,:,bR[a],xl] = f[:,:,biR[a],xl]
              g[:,:,bR[a],xl] = g[:,:,biR[a],xl]
              if rheology == "archean2":
                 h[:,:,bR[a],xl] = h[:,:,biR[a],xl]

          for a in range(naR,nbR):   
              tmp              = 0.5*(f[:,:,bR[a],xl] + f[:,:,biR[a],xl])
              f[:,:,bR [a],xl] = tmp
              f[:,:,biR[a],xl] = tmp
              tmp              = 0.5*(g[:,:,bR[a],xl] + g[:,:,biR[a],xl])
              g[:,:,bR [a],xl] = tmp
              g[:,:,biR[a],xl] = tmp
              if rheology == "archean2":
                 tmp              = 0.5*(h[:,:,bR[a],xl] + h[:,:,biR[a],xl])
                 h[:,:,bR [a],xl] = tmp
                 h[:,:,biR[a],xl] = tmp


    if boundary_r:
       if RKTLBM:  
          for a in range(naL):
              f_r[:,:,bL[a],xr] = f_r[:,:,biL[a],xr]
              f_b[:,:,bL[a],xr] = f_b[:,:,biL[a],xr]
              g[:,:,bL[a],xr]   = g  [:,:,biL[a],xr]
              if rheology == "archean2":
                 h[:,:,bL[a],xr] = h[:,:,biL[a],xr]

          for a in range(naL,nbL):
              tmp              = 0.5*(f_r[:,:,bL[a],xr] + f_r[:,:,biL[a],xr])
              f_r[:,:,bL [a],xr] = tmp
              f_r[:,:,biL[a],xr] = tmp
              tmp              = 0.5*(f_b[:,:,bL[a],xr] + f_b[:,:,biL[a],xr])
              f_b[:,:,bL [a],xr] = tmp
              f_b[:,:,biL[a],xr] = tmp
              tmp              = 0.5*(g[:,:,bL[a],xr] + g[:,:,biL[a],xr])
              g[:,:,bL [a],xr] = tmp
              g[:,:,biL[a],xr] = tmp
              if rheology == "archean2":
                 tmp              = 0.5*(h[:,:,bL[a],xr] + h[:,:,biL[a],xr])
                 h[:,:,bL [a],xr] = tmp
                 h[:,:,biL[a],xr] = tmp

       else:
          for a in range(naL):
              f[:,:,bL[a],xr] = f[:,:,biL[a],xr]
              g[:,:,bL[a],xr] = g[:,:,biL[a],xr]
              if rheology == "archean2":
                 h[:,:,bL[a],xr] = h[:,:,biL[a],xr]

          for a in range(naL,nbL):
              tmp              = 0.5*(f[:,:,bL[a],xr] + f[:,:,biL[a],xr])
              f[:,:,bL [a],xr] = tmp
              f[:,:,biL[a],xr] = tmp
              tmp              = 0.5*(g[:,:,bL[a],xr] + g[:,:,biL[a],xr])
              g[:,:,bL [a],xr] = tmp
              g[:,:,biL[a],xr] = tmp
              if rheology == "archean2":
                 tmp              = 0.5*(h[:,:,bL[a],xr] + h[:,:,biL[a],xr])
                 h[:,:,bL [a],xr] = tmp
                 h[:,:,biL[a],xr] = tmp



#   Upper and lower BC

    if boundary_u:
        if RKTLBM:
           for a in range(naD):
               f_r[:,zu,bD[a],:] = f_r[:,zu,biD[a],:]
               f_b[:,zu,bD[a],:] = f_b[:,zu,biD[a],:]
               g[:,zu,bD[a],:]   = g  [:,zu,biD[a],:]
               if rheology == "archean2":
                  h[:,zu,bD[a],:] = h[:,zu,biD[a],:]
           for a in range(naD,nbD):
               tmp              = 0.5*(f_r[:,zu,bD[a],:] + f_r[:,zu,biD[a],:])
               f_r[:,zu,bD [a],:] = tmp
               f_r[:,zu,biD[a],:] = tmp
               tmp              = 0.5*(f_b[:,zu,bD[a],:] + f_b[:,zu,biD[a],:])
               f_b[:,zu,bD [a],:] = tmp
               f_b[:,zu,biD[a],:] = tmp
               tmp              = 0.5*(g[:,zu,bD[a],:] + g[:,zu,biD[a],:])
               g[:,zu,bD [a],:] = tmp
               g[:,zu,biD[a],:] = tmp
               if rheology == "archean2":
                  tmp              = 0.5*(h[:,zu,bD[a],:] + h[:,zu,biD[a],:])
                  h[:,zu,bD [a],:] = tmp
                  h[:,zu,biD[a],:] = tmp

        else:
           for a in range(naD):
               f[:,zu,bD[a],:] = f[:,zu,biD[a],:]
               g[:,zu,bD[a],:] = g[:,zu,biD[a],:]
               if rheology == "archean2":
                  h[:,zu,bD[a],:] = h[:,zu,biD[a],:]
           for a in range(naD,nbD):
               tmp              = 0.5*(f[:,zu,bD[a],:] + f[:,zu,biD[a],:])
               f[:,zu,bD [a],:] = tmp
               f[:,zu,biD[a],:] = tmp
               tmp              = 0.5*(g[:,zu,bD[a],:] + g[:,zu,biD[a],:])
               g[:,zu,bD [a],:] = tmp
               g[:,zu,biD[a],:] = tmp
               if rheology == "archean2":
                  tmp              = 0.5*(h[:,zu,bD[a],:] + h[:,zu,biD[a],:])
                  h[:,zu,bD [a],:] = tmp
                  h[:,zu,biD[a],:] = tmp

    if boundary_d:
       if RKTLBM:
          for a in range(naU):
              f_r[:,zd,bU[a],:] = f_r[:,zd,biU[a],:]
              f_b[:,zd,bU[a],:] = f_b[:,zd,biU[a],:]
              g[:,zd,bU[a],:]   = g  [:,zd,biU[a],:]
              if rheology == "archean2":
                 h[:,zd,bU[a],:] = h[:,zd,biU[a],:]
          for a in range(naU,nbU):
              tmp              = 0.5*(f_r[:,zd,bU[a],:] + f_r[:,zd,biU[a],:])
              f_r[:,zd,bU [a],:] = tmp
              f_r[:,zd,biU[a],:] = tmp
              tmp              = 0.5*(f_b[:,zd,bU[a],:] + f_b[:,zd,biU[a],:])
              f_b[:,zd,bU [a],:] = tmp
              f_b[:,zd,biU[a],:] = tmp
              tmp              = 0.5*(g[:,zd,bU[a],:] + g[:,zd,biU[a],:])
              g[:,zd,bU [a],:] = tmp
              g[:,zd,biU[a],:] = tmp
              if rheology == "archean2":
                 tmp              = 0.5*(h[:,zd,bU[a],:] + h[:,zd,biU[a],:])
                 h[:,zd,bU [a],:] = tmp
                 h[:,zd,biU[a],:] = tmp

       else:  
          for a in range(naU):
              f[:,zd,bU[a],:] = f[:,zd,biU[a],:]
              g[:,zd,bU[a],:] = g[:,zd,biU[a],:]
              if rheology == "archean2":
                 h[:,zd,bU[a],:] = h[:,zd,biU[a],:]
          for a in range(naU,nbU):
              tmp              = 0.5*(f[:,zd,bU[a],:] + f[:,zd,biU[a],:])
              f[:,zd,bU [a],:] = tmp
              f[:,zd,biU[a],:] = tmp
              tmp              = 0.5*(g[:,zd,bU[a],:] + g[:,zd,biU[a],:])
              g[:,zd,bU [a],:] = tmp
              g[:,zd,biU[a],:] = tmp
              if rheology == "archean2":
                 tmp              = 0.5*(h[:,zd,bU[a],:] + h[:,zd,biU[a],:])
                 h[:,zd,bU [a],:] = tmp
                 h[:,zd,biU[a],:] = tmp


#   Forward and backward BC

    if D == 3:
      if boundary_b:
         if RKTLBM:
            for a in range(naF):
                f_r[yb,:,bF[a],:] = f_r[yb,:,biF[a],:]
                f_b[yb,:,bF[a],:] = f_b[yb,:,biF[a],:]
                g[yb,:,bF[a],:]   = g  [yb,:,biF[a],:]
                if rheology == "archean2":
                   h[yb,:,bF[a],:] = h[yb,:,biF[a],:]
            for a in range(naF,nbF):
                tmp              = 0.5*(f_r[yb,:,bF[a],:] + f_r[yb,:,biF[a],:])
                f_r[yb,:,bF [a],:] = tmp
                f_r[yb,:,biF[a],:] = tmp
                tmp              = 0.5*(f_b[yb,:,bF[a],:] + f_b[yb,:,biF[a],:])
                f_b[yb,:,bF [a],:] = tmp
                f_b[yb,:,biF[a],:] = tmp
                tmp              = 0.5*(g[yb,:,bF[a],:] + g[yb,:,biF[a],:])
                g[yb,:,bF [a],:] = tmp
                g[yb,:,biF[a],:] = tmp
                if rheology == "archean2":
                   tmp              = 0.5*(h[yb,:,bF[a],:] + h[yb,:,biF[a],:])
                   h[yb,:,bF [a],:] = tmp
                   h[yb,:,biF[a],:] = tmp
         
         else:
            for a in range(naF):
                f[yb,:,bF[a],:] = f[yb,:,biF[a],:]
                g[yb,:,bF[a],:] = g[yb,:,biF[a],:]
                if rheology == "archean2":
                   h[yb,:,bF[a],:] = h[yb,:,biF[a],:]
            for a in range(naF,nbF):
                tmp              = 0.5*(f[yb,:,bF[a],:] + f[yb,:,biF[a],:])
                f[yb,:,bF [a],:] = tmp
                f[yb,:,biF[a],:] = tmp
                tmp              = 0.5*(g[yb,:,bF[a],:] + g[yb,:,biF[a],:])
                g[yb,:,bF [a],:] = tmp
                g[yb,:,biF[a],:] = tmp
                if rheology == "archean2":
                   tmp              = 0.5*(h[yb,:,bF[a],:] + h[yb,:,biF[a],:])
                   h[yb,:,bF [a],:] = tmp
                   h[yb,:,biF[a],:] = tmp

      if boundary_f:
         if RKTLBM:
            for a in range(naB):
                f_r[yf,:,bB[a],:] = f_r[yf,:,biB[a],:]
                f_b[yf,:,bB[a],:] = f_b[yf,:,biB[a],:]
                g[yf,:,bB[a],:]   = g  [yf,:,biB[a],:]
                if rheology == "archean2":
                   h[yf,:,bB[a],:] = h[yf,:,biB[a],:]
            for a in range(naB,nbB):
                tmp              = 0.5*(f_r[yb,:,bB[a],:] + f_r[yb,:,biB[a],:])
                f_r[yb,:,bB [a],:] = tmp
                f_r[yb,:,biB[a],:] = tmp
                tmp              = 0.5*(f_b[yb,:,bB[a],:] + f_b[yb,:,biB[a],:])
                f_b[yb,:,bB [a],:] = tmp
                f_b[yb,:,biB[a],:] = tmp
                tmp              = 0.5*(g[yb,:,bB[a],:] + g[yb,:,biB[a],:])
                g[yb,:,bB [a],:] = tmp
                g[yb,:,biB[a],:] = tmp
                if rheology == "archean2":
                   tmp              = 0.5*(h[yb,:,bB[a],:] + h[yb,:,biB[a],:])
                   h[yb,:,bB [a],:] = tmp
                   h[yb,:,biB[a],:] = tmp

         else:
            for a in range(naB):
                f[yf,:,bB[a],:] = f[yf,:,biB[a],:]
                g[yf,:,bB[a],:] = g[yf,:,biB[a],:]
                if rheology == "archean2":
                   h[yf,:,bB[a],:] = h[yf,:,biB[a],:]
            for a in range(naB,nbB):
                tmp              = 0.5*(f[yb,:,bB[a],:] + f[yb,:,biB[a],:])
                f[yb,:,bB [a],:] = tmp
                f[yb,:,biB[a],:] = tmp
                tmp              = 0.5*(g[yb,:,bB[a],:] + g[yb,:,biB[a],:])
                g[yb,:,bB [a],:] = tmp
                g[yb,:,biB[a],:] = tmp
                if rheology == "archean2":
                   tmp              = 0.5*(h[yb,:,bB[a],:] + h[yb,:,biB[a],:])
                   h[yb,:,bB [a],:] = tmp
                   h[yb,:,biB[a],:] = tmp


#===============================================================================
# Get the edges to the regions in each node
#===============================================================================

#-------------------------------------------------------------------------------
# Initialization get_edges functions that minimize overhead and are fast
# by using persistent send/recv requests
#-------------------------------------------------------------------------------

def initialize_get_edges():

    global i_F, i_G, i_H, i_F_r, i_F_b, i_HD, i_HB
    global i_rho_r, i_rho_b, i_D_ed, i_D_eb
    global    aL,    aR,    aU,    aD,    aB,    aF
    global   naL,   naR,   naU,   naD,   naB,   naF
    global edgeL, edgeR, edgeU, edgeD, edgeB, edgeF
    global edge_rho_L, edge_rho_R, edge_rho_U, edge_rho_D, edge_rho_B, edge_rho_F
    global   fgL    ,   fgR    ,   fgU    ,   fgD    ,   fgB    ,   fgF
    global   rhoL   ,   rhoR   ,   rhoU   ,   rhoD   ,   rhoB   ,   rhoF
    global R, Rx, Ry, Rz, rr, rl, rf, rb, rd, ru
    global zcol_comm
    global X, Y, Z

#   Faces

    global reqSendR, reqSendL, reqRecvR, reqRecvL
    global reqSendD, reqSendU, reqRecvD, reqRecvU
    global reqSendF, reqSendB, reqRecvF, reqRecvB
    global reqSendrhoR, reqSendrhoL, reqRecvrhoR, reqRecvrhoL
    global reqSendrhoD, reqSendrhoU, reqRecvrhoD, reqRecvrhoU
    global reqSendrhoF, reqSendrhoB, reqRecvrhoF, reqRecvrhoB

#   Corners of cuboid

    global    edgeRD,    edgeRU,    edgeLD,    edgeLU
    global      fgRD,      fgRU,      fgLD,      fgLU
    global      naRD,      naRU,      naLD,      naLU
    global rrd,  aRD, rru,  aRU, rld,  aLD, rlu,  aLU
    global reqSendRD, reqSendRU, reqSendLD, reqSendLU
    global reqRecvRD, reqRecvRU, reqRecvLD, reqRecvLU

    global    edgeFD,    edgeFU,    edgeBD,    edgeBU
    global      fgFD,      fgFU,      fgBD,      fgBU
    global      naFD,      naFU,      naBD,      naBU
    global rfd,  aFD, rfu,  aFU, rbd,  aBD, rbu,  aBU
    global reqSendFD, reqSendFU, reqSendBD, reqSendBU
    global reqRecvFD, reqRecvFU, reqRecvBD, reqRecvBU

    global    edgeRF,    edgeRB,    edgeLF,    edgeLB
    global      fgRF,      fgRB,      fgLF,      fgLB
    global      naRF,      naRB,      naLF,      naLB
    global rrf,  aRF, rrb,  aRB, rlf,  aLF, rlb,  aLB
    global reqSendRF, reqSendRB, reqSendLF, reqSendLB
    global reqRecvRF, reqRecvRB, reqRecvLF, reqRecvLB

    global xedge, yedge, zedge, nx_i, ny_i, nz_i, mxi, myi, mzi, mx_i, my_i, mz_i

#   Set the X, Y and Z indices for 2D and 3D problems (Z index is direction of gravity)

    if D == 2:
        X = 0
        Y = 1
        Z = 1
    else:
        X = 0
        Y = 1
        Z = 2

    xedge = 1 ; yedge=1 ; zedge = 1
    if nrx == 1: xedge = 0
    if nry == 1: yedge = 0
    if nrz == 1: zedge = 0

    nx_i = nxi+2*xedge
    ny_i = nyi+2*yedge
    nz_i = nzi+2*zedge

    R  = size
    Rx = size
    Ry = size
    Rz = size
    for rx in range(nrx):
        for ry in range(nry):
            for rz in range(nrz):
                r = rx + rz*nrx + ry*nrx*nrz
                if rank == r:
                    R    = r
                    Rx   = rx
                    Ry   = ry
                    Rz   = rz

                    mxi  = mx[rx]
                    myi  = my[ry]
                    mzi  = mz[rz]
                    mx_i = mxi + 2*xedge
                    my_i = myi + 2*yedge
                    mz_i = mzi + 2*zedge

                    rx_r = (rx+1    )%nrx
                    rx_l = (rx-1+nrx)%nrx
                    rr   = rx_r+rz*nrx+ry*nrx*nrz
                    rl   = rx_l+rz*nrx+ry*nrx*nrz

                    rz_d = (rz+1    )%nrz
                    rz_u = (rz-1+nrz)%nrz
                    rd   = rx+rz_d*nrx+ry*nrx*nrz
                    ru   = rx+rz_u*nrx+ry*nrx*nrz

                    ry_f = (ry+1    )%nry
                    ry_b = (ry-1+nry)%nry
                    rf   = rx+rz*nrx+ry_f*nrx*nrz
                    rb   = rx+rz*nrx+ry_b*nrx*nrz

                    rlu  = rx_l +rz_u*nrx + ry*nrx*nrz
                    rld  = rx_l +rz_d*nrx + ry*nrx*nrz
                    rru  = rx_r +rz_u*nrx + ry*nrx*nrz
                    rrd  = rx_r +rz_d*nrx + ry*nrx*nrz

                    rlb  = rx_l +rz*nrx + ry_b*nrx*nrz
                    rlf  = rx_l +rz*nrx + ry_f*nrx*nrz
                    rrb  = rx_r +rz*nrx + ry_b*nrx*nrz
                    rrf  = rx_r +rz*nrx + ry_f*nrx*nrz

                    rbu  = rx +rz_u*nrx + ry_b*nrx*nrz
                    rbd  = rx +rz_d*nrx + ry_b*nrx*nrz
                    rfu  = rx +rz_u*nrx + ry_f*nrx*nrz
                    rfd  = rx +rz_d*nrx + ry_f*nrx*nrz

    # Z-column communicator: groups all ranks sharing the same (Rx, Ry) but
    # different Rz.  Used for efficient melt-surface transport (avoids full
    # field gather/scatter).  Within each column comm, Rz==0 is root (rank 0).
    zcol_color = (Rx + Ry * nrx) if rank < nr else MPI.UNDEFINED
    zcol_comm  = comm.Split(zcol_color, Rz if rank < nr else 0)

    if RKTLBM:
       i_F_r    = 0
       i_rho_r  = 0
       i_F_b    = 1
       i_rho_b  = 1
       i_D_ed   = 2   # advected depletion ghosts ride the rho exchange
       i_D_eb   = 3
       i_G      = 2
       i_H      = 3
    else: 
       i_F      = 0
       i_G      = 1
       i_H      = 2

    aL   = np.zeros(na,dtype=np_int)
    aR   = np.zeros(na,dtype=np_int)
    aU   = np.zeros(na,dtype=np_int)
    aD   = np.zeros(na,dtype=np_int)
    aB   = np.zeros(na,dtype=np_int)
    aF   = np.zeros(na,dtype=np_int)
    naL  = np_int(0) ; naR  = np_int(0)
    naU  = np_int(0) ; naD  = np_int(0)
    naB  = np_int(0) ; naF  = np_int(0)
    naLU = 0 ; naLD = 0 ; naRU = 0 ; naRD = 0
    naBU = 0 ; naBD = 0 ; naFU = 0 ; naFD = 0
    naLB = 0 ; naLF = 0 ; naRB = 0 ; naRF = 0
    for a in range(1,na):
        if   c[a,X] < 0: # Moving left
            aL[naL] = a
            naL    += 1
        elif c[a,X] > 0: # Moving right
            aR[naR] = a
            naR    += 1
        if D == 3:
            if   c[a,Y] < 0: # Moving backward
                aB[naB] = a
                naB    += 1
            elif c[a,Y] > 0: # Moving forward
                aF[naF] = a
                naF    += 1
        if   c[a,Z] < 0: # Moving up
            aU[naU] = a
            naU    += 1
        elif c[a,Z] > 0: # Moving down
            aD[naD] = a
            naD    += 1

        # Corners

        if xedge == 1 and zedge == 1:
            if   c[a,X] ==  1 and c[a,Z] ==  1:
                naRD = 1 ; aRD = a
            if   c[a,X] ==  1 and c[a,Z] == -1:
                naRU = 1 ; aRU = a
            if   c[a,X] == -1 and c[a,Z] ==  1:
                naLD = 1 ; aLD = a
            if   c[a,X] == -1 and c[a,Z] == -1:
                naLU = 1 ; aLU = a

        if yedge == 1 and zedge == 1:
            if   c[a,Y] ==  1 and c[a,Z] ==  1:
                naFD = 1 ; aFD = a
            if   c[a,Y] ==  1 and c[a,Z] == -1:
                naFU = 1 ; aFU = a
            if   c[a,Y] == -1 and c[a,Z] ==  1:
                naBD = 1 ; aBD = a
            if   c[a,Y] == -1 and c[a,Z] == -1:
                naBU = 1 ; aBU = a

        if xedge == 1 and yedge == 1:
            if   c[a,X] ==  1 and c[a,Y] ==  1:
                naRF = 1 ; aRF = a
            if   c[a,X] ==  1 and c[a,Y] == -1:
                naRB = 1 ; aRB = a
            if   c[a,X] == -1 and c[a,Y] ==  1:
                naLF = 1 ; aLF = a
            if   c[a,X] == -1 and c[a,Y] == -1:
                naLB = 1 ; aLB = a

#   Setup buffers for persistent sends and receives

    if rank < nr:

        if not (overlap and code == "c"): # Case of no overlapping
            if RKTLBM :
               # hd/hb ride the packed exchange for archean2 (slots 4,5)
               nfg_rk      = 4   # hd/hb retired from the packed exchange (TVD scalars now)
               fgL         = np.zeros((nfg_rk,naL,my_i,mz_i),dtype=np_edge)
               fgR         = np.zeros((nfg_rk,naL,my_i,mz_i),dtype=np_edge)
               fgU         = np.zeros((nfg_rk,naU,my_i,mx_i),dtype=np_edge)
               fgD         = np.zeros((nfg_rk,naU,my_i,mx_i),dtype=np_edge)
               fgB         = np.zeros((nfg_rk,naB,mz_i,mx_i),dtype=np_edge)
               fgF         = np.zeros((nfg_rk,naF,mz_i,mx_i),dtype=np_edge)

               edgeL       = np.zeros((nfg_rk,naL,my_i,mz_i),dtype=np_edge)
               edgeR       = np.zeros((nfg_rk,naL,my_i,mz_i),dtype=np_edge)
               edgeU       = np.zeros((nfg_rk,naU,my_i,mx_i),dtype=np_edge)
               edgeD       = np.zeros((nfg_rk,naU,my_i,mx_i),dtype=np_edge)
               edgeB       = np.zeros((nfg_rk,naB,mz_i,mx_i),dtype=np_edge)
               edgeF       = np.zeros((nfg_rk,naF,mz_i,mx_i),dtype=np_edge)

               rhoL        = np.zeros((4,my_i,mz_i),dtype=np_edge)
               rhoR        = np.zeros((4,my_i,mz_i),dtype=np_edge)
               rhoU        = np.zeros((2,my_i,mx_i),dtype=np_edge)
               rhoD        = np.zeros((2,my_i,mx_i),dtype=np_edge)
               rhoB        = np.zeros((2,mz_i,mx_i),dtype=np_edge)
               rhoF        = np.zeros((2,mz_i,mx_i),dtype=np_edge)

               edge_rho_L  = np.zeros((4,my_i,mz_i),dtype=np_edge)
               edge_rho_R  = np.zeros((4,my_i,mz_i),dtype=np_edge)
               edge_rho_U  = np.zeros((2,my_i,mx_i),dtype=np_edge)
               edge_rho_D  = np.zeros((2,my_i,mx_i),dtype=np_edge)
               edge_rho_B  = np.zeros((2,mz_i,mx_i),dtype=np_edge)
               edge_rho_F  = np.zeros((2,mz_i,mx_i),dtype=np_edge)

#           Faces including edges
            else:
               fgL         = np.zeros((3,naL,my_i,mz_i),dtype=np_edge)
               fgR         = np.zeros((3,naL,my_i,mz_i),dtype=np_edge)
               fgU         = np.zeros((3,naU,my_i,mx_i),dtype=np_edge)
               fgD         = np.zeros((3,naU,my_i,mx_i),dtype=np_edge)
               fgB         = np.zeros((3,naB,mz_i,mx_i),dtype=np_edge)
               fgF         = np.zeros((3,naF,mz_i,mx_i),dtype=np_edge)

               edgeL       = np.zeros((3,naL,my_i,mz_i),dtype=np_edge)
               edgeR       = np.zeros((3,naL,my_i,mz_i),dtype=np_edge)
               edgeU       = np.zeros((3,naU,my_i,mx_i),dtype=np_edge)
               edgeD       = np.zeros((3,naU,my_i,mx_i),dtype=np_edge)
               edgeB       = np.zeros((3,naB,mz_i,mx_i),dtype=np_edge)
               edgeF       = np.zeros((3,naF,mz_i,mx_i),dtype=np_edge)
            
               rhoL        = np.zeros((4,my_i,mz_i),dtype=np_edge)
               rhoR        = np.zeros((4,my_i,mz_i),dtype=np_edge)
               rhoU        = np.zeros((2,my_i,mx_i),dtype=np_edge)
               rhoD        = np.zeros((2,my_i,mx_i),dtype=np_edge)
               rhoB        = np.zeros((2,mz_i,mx_i),dtype=np_edge)
               rhoF        = np.zeros((2,mz_i,mx_i),dtype=np_edge)
            
               edge_rho_L  = np.zeros((4,my_i,mz_i),dtype=np_edge)
               edge_rho_R  = np.zeros((4,my_i,mz_i),dtype=np_edge)
               edge_rho_U  = np.zeros((2,my_i,mx_i),dtype=np_edge)
               edge_rho_D  = np.zeros((2,my_i,mx_i),dtype=np_edge)
               edge_rho_B  = np.zeros((2,mz_i,mx_i),dtype=np_edge)
               edge_rho_F  = np.zeros((2,mz_i,mx_i),dtype=np_edge)

        else:           # Case of overlapping

#           Faces excluding edges

            fgL     = np.zeros((3,naL,myi,mzi),dtype=np_edge)
            fgR     = np.zeros((3,naL,myi,mzi),dtype=np_edge)
            fgU     = np.zeros((3,naU,myi,mxi),dtype=np_edge)
            fgD     = np.zeros((3,naU,myi,mxi),dtype=np_edge)
            fgB     = np.zeros((3,naB,mzi,mxi),dtype=np_edge)
            fgF     = np.zeros((3,naF,mzi,mxi),dtype=np_edge)

            edgeL   = np.zeros((3,naL,myi,mzi),dtype=np_edge)
            edgeR   = np.zeros((3,naL,myi,mzi),dtype=np_edge)
            edgeU   = np.zeros((3,naU,myi,mxi),dtype=np_edge)
            edgeD   = np.zeros((3,naU,myi,mxi),dtype=np_edge)
            edgeB   = np.zeros((3,naB,mzi,mxi),dtype=np_edge)
            edgeF   = np.zeros((3,naF,mzi,mxi),dtype=np_edge)

            # Edges normal to y-plane (ie. xz corners)

            if naLU == 1:

                fgLU   = np.zeros((3,myi),dtype=np_edge)
                fgLD   = np.zeros((3,myi),dtype=np_edge)
                fgRU   = np.zeros((3,myi),dtype=np_edge)
                fgRD   = np.zeros((3,myi),dtype=np_edge)

                edgeLU = np.zeros((3,myi),dtype=np_edge)
                edgeLD = np.zeros((3,myi),dtype=np_edge)
                edgeRU = np.zeros((3,myi),dtype=np_edge)
                edgeRD = np.zeros((3,myi),dtype=np_edge)

            # Edges normal to x-plane (ie. yz corners)

            if naBU == 1:
                fgBU   = np.zeros((3,mxi),dtype=np_edge)
                fgBD   = np.zeros((3,mxi),dtype=np_edge)
                fgFU   = np.zeros((3,mxi),dtype=np_edge)
                fgFD   = np.zeros((3,mxi),dtype=np_edge)

                edgeBU = np.zeros((3,mxi),dtype=np_edge)
                edgeBD = np.zeros((3,mxi),dtype=np_edge)
                edgeFU = np.zeros((3,mxi),dtype=np_edge)
                edgeFD = np.zeros((3,mxi),dtype=np_edge)

            # Edges normal to z-plane (ie. xy corners)

            if naLB == 1:
                fgLB   = np.zeros((3,mzi),dtype=np_edge)
                fgLF   = np.zeros((3,mzi),dtype=np_edge)
                fgRB   = np.zeros((3,mzi),dtype=np_edge)
                fgRF   = np.zeros((3,mzi),dtype=np_edge)

                edgeLB = np.zeros((3,mzi),dtype=np_edge)
                edgeLF = np.zeros((3,mzi),dtype=np_edge)
                edgeRB = np.zeros((3,mzi),dtype=np_edge)
                edgeRF = np.zeros((3,mzi),dtype=np_edge)

#       Persistent send and receive requests for faces

        if nrx > 1: # X persistent send and receive requests for faces
            reqSendR = comm.Send_init([fgR  ,MPI_EDGE] ,dest  =rl,tag=1)
            reqSendL = comm.Send_init([fgL  ,MPI_EDGE] ,dest  =rr,tag=2)
            reqRecvR = comm.Recv_init([edgeR,MPI_EDGE] ,source=rr,tag=1)
            reqRecvL = comm.Recv_init([edgeL,MPI_EDGE] ,source=rl,tag=2)
            if RKTLBM :
               reqSendrhoR = comm.Send_init([rhoR  ,MPI_EDGE] ,dest  =rl,tag=7)
               reqSendrhoL = comm.Send_init([rhoL  ,MPI_EDGE] ,dest  =rr,tag=8)
               reqRecvrhoR = comm.Recv_init([edge_rho_R,MPI_EDGE] ,source=rr,tag=7)
               reqRecvrhoL = comm.Recv_init([edge_rho_L,MPI_EDGE] ,source=rl,tag=8)
        if nrz > 1: # Z persistent send and receive requests for faces
            reqSendD = comm.Send_init([fgD  ,MPI_EDGE] ,dest  =ru,tag=3)
            reqSendU = comm.Send_init([fgU  ,MPI_EDGE] ,dest  =rd,tag=4)
            reqRecvD = comm.Recv_init([edgeD,MPI_EDGE] ,source=rd,tag=3)
            reqRecvU = comm.Recv_init([edgeU,MPI_EDGE] ,source=ru,tag=4)
            if RKTLBM:
               reqSendrhoD = comm.Send_init([rhoD  ,MPI_EDGE] ,dest  =ru,tag=9)
               reqSendrhoU = comm.Send_init([rhoU  ,MPI_EDGE] ,dest  =rd,tag=10)
               reqRecvrhoD = comm.Recv_init([edge_rho_D,MPI_EDGE] ,source=rd,tag=9)
               reqRecvrhoU = comm.Recv_init([edge_rho_U,MPI_EDGE] ,source=ru,tag=10)
        if nry > 1: # Y persistent send and receive requests for faces
            reqSendF = comm.Send_init([fgF  ,MPI_EDGE] ,dest  =rb,tag=5)
            reqSendB = comm.Send_init([fgB  ,MPI_EDGE] ,dest  =rf,tag=6)
            reqRecvF = comm.Recv_init([edgeF,MPI_EDGE] ,source=rf,tag=5)
            reqRecvB = comm.Recv_init([edgeB,MPI_EDGE] ,source=rb,tag=6)
            if RKTLBM:
               reqSendrhoF = comm.Send_init([rhoF  ,MPI_EDGE] ,dest  =rb,tag=11)
               reqSendrhoB = comm.Send_init([rhoB  ,MPI_EDGE] ,dest  =rf,tag=12)
               reqRecvrhoF = comm.Recv_init([edge_rho_F,MPI_EDGE] ,source=rf,tag=11)
               reqRecvrhoB = comm.Recv_init([edge_rho_B,MPI_EDGE] ,source=rb,tag=12)

#       Persistent send and receove requests for face corners needed for overlapping

        if overlap and code == "c":

            if naRD == 1: # Corners of XZ faces
                reqSendRD = comm.Send_init([fgRD  ,MPI_EDGE] ,dest  =rlu,tag=11)
                reqSendRU = comm.Send_init([fgRU  ,MPI_EDGE] ,dest  =rld,tag=12)
                reqSendLD = comm.Send_init([fgLD  ,MPI_EDGE] ,dest  =rru,tag=13)
                reqSendLU = comm.Send_init([fgLU  ,MPI_EDGE] ,dest  =rrd,tag=14)

                reqRecvRD = comm.Recv_init([edgeRD,MPI_EDGE] ,source=rrd,tag=11)
                reqRecvRU = comm.Recv_init([edgeRU,MPI_EDGE] ,source=rru,tag=12)
                reqRecvLD = comm.Recv_init([edgeLD,MPI_EDGE] ,source=rld,tag=13)
                reqRecvLU = comm.Recv_init([edgeLU,MPI_EDGE] ,source=rlu,tag=14)

            if naFD == 1: # Corners of YZ faces
                reqSendFD = comm.Send_init([fgFD  ,MPI_EDGE] ,dest  =rbu,tag=21)
                reqSendFU = comm.Send_init([fgFU  ,MPI_EDGE] ,dest  =rbd,tag=22)
                reqSendBD = comm.Send_init([fgBD  ,MPI_EDGE] ,dest  =rfu,tag=23)
                reqSendBU = comm.Send_init([fgBU  ,MPI_EDGE] ,dest  =rfd,tag=24)

                reqRecvFD = comm.Recv_init([edgeFD,MPI_EDGE] ,source=rfd,tag=21)
                reqRecvFU = comm.Recv_init([edgeFU,MPI_EDGE] ,source=rfu,tag=22)
                reqRecvBD = comm.Recv_init([edgeBD,MPI_EDGE] ,source=rbd,tag=23)
                reqRecvBU = comm.Recv_init([edgeBU,MPI_EDGE] ,source=rbu,tag=24)

            if naRF == 1: # Corners of XY faces
                reqSendRF = comm.Send_init([fgRF  ,MPI_EDGE] ,dest  =rlb,tag=31)
                reqSendRB = comm.Send_init([fgRB  ,MPI_EDGE] ,dest  =rlf,tag=32)
                reqSendLF = comm.Send_init([fgLF  ,MPI_EDGE] ,dest  =rrb,tag=33)
                reqSendLB = comm.Send_init([fgLB  ,MPI_EDGE] ,dest  =rrf,tag=34)

                reqRecvRF = comm.Recv_init([edgeRF,MPI_EDGE] ,source=rrf,tag=31)
                reqRecvRB = comm.Recv_init([edgeRB,MPI_EDGE] ,source=rrb,tag=32)
                reqRecvLF = comm.Recv_init([edgeLF,MPI_EDGE] ,source=rlf,tag=33)
                reqRecvLB = comm.Recv_init([edgeLB,MPI_EDGE] ,source=rlb,tag=34)

#-------------------------------------------------------------------------------
# Get of the edges for the case overlap = False
#-------------------------------------------------------------------------------

def get_edges():

    if nrx > 1 and rank < nr:
        get_edgesX()   # Get X edges
    if nrz > 1 and rank < nr:
        get_edgesZ()   # Get Z edges
    if nry > 1 and rank < nr:
        get_edgesY()   # Get Y edges

#-------------------------------------------------------------------------------
# Kernel to get the X edges of a domain for overlap = False
#-------------------------------------------------------------------------------

def get_edgesX():

    global i_F, i_G, i_H, i_F_r, i_F_b, i_HD, i_HB
    global i_rho_r, i_rho_b, i_D_ed, i_D_eb
    global    aL,    aR,    aU,    aD,    aB,    aF
    global   naL,   naR,   naU,   naD,   naB,   naF
    global edgeL    , edgeR    , edgeU    , edgeD    , edgeB    , edgeF
    global edge_rho_L    , edge_rho_R    , edge_rho_U    , edge_rho_D    , edge_rho_B    , edge_rho_F
    global   fgL    ,   fgR    ,   fgU    ,   fgD    ,   fgB    ,   fgF
    global   rhoL    ,  rhoR    ,  rhoU    ,  rhoD    ,  rhoB    ,  rhoF
    global R, Rx, Ry, Rz, rr, rl, rf, rb, rd, ru
    global reqSendR, reqSendL, reqRecvR, reqRecvL
    global reqSendD, reqSendU, reqRecvD, reqRecvU
    global reqSendF, reqSendB, reqRecvF, reqRecvB

    if RKTLBM:
       for a in range(naL):
           fgR[i_F_r,a] = f_r [0:my_i,0:mz_i,aL[a],1     ]
           fgR[i_F_b,a] = f_b [0:my_i,0:mz_i,aL[a],1     ]
           fgR[i_G,a]   = g   [0:my_i,0:mz_i,aL[a],1     ]
           if rheology == "archean2":
              fgR[i_H,a] = h  [0:my_i,0:mz_i,aL[a],1     ]
    else:
       for a in range(naL):
           fgR[i_F,a] = f[0:my_i,0:mz_i,aL[a],1     ]
           fgR[i_G,a] = g[0:my_i,0:mz_i,aL[a],1     ]
           if rheology == "archean2":
              fgR[i_H,a] = h[0:my_i,0:mz_i,aL[a],1     ]

    if RKTLBM:
       for a in range(naR):
           fgL[i_F_r,a] = f_r [0:my_i,0:mz_i,aR[a],mx[Rx]]
           fgL[i_F_b,a] = f_b [0:my_i,0:mz_i,aR[a],mx[Rx]]
           fgL[i_G,a]   = g   [0:my_i,0:mz_i,aR[a],mx[Rx]]
           if rheology == "archean2":
              fgL[i_H,a] = h  [0:my_i,0:mz_i,aR[a],mx[Rx]]
           
    else:       
       for a in range(naR):
           fgL[i_F,a] = f[0:my_i,0:mz_i,aR[a],mx[Rx]]
           fgL[i_G,a] = g[0:my_i,0:mz_i,aR[a],mx[Rx]]
           if rheology == "archean2":
              fgL[i_H,a] = h[0:my_i,0:mz_i,aR[a],mx[Rx]]

    MPI.Prequest.Startall([reqSendR, reqRecvR, reqSendL, reqRecvL])
    MPI.Request.Waitall  ([reqSendR, reqRecvR, reqSendL, reqRecvL])

    if RKTLBM: 
       for a in range(naL):
           f_r[0:my_i,0:mz_i,aL[a],mx[Rx]+1] = edgeR[i_F_r,a]
           f_b[0:my_i,0:mz_i,aL[a],mx[Rx]+1] = edgeR[i_F_b,a]
           g  [0:my_i,0:mz_i,aL[a],mx[Rx]+1] = edgeR[i_G,a]
           if rheology == "archean2":
              h[0:my_i,0:mz_i,aL[a],mx[Rx]+1] = edgeR[i_H,a]
    else:
       for a in range(naL):
           f[0:my_i,0:mz_i,aL[a],mx[Rx]+1] = edgeR[i_F,a]
           g[0:my_i,0:mz_i,aL[a],mx[Rx]+1] = edgeR[i_G,a]
           if rheology == "archean2":
              h[0:my_i,0:mz_i,aL[a],mx[Rx]+1] = edgeR[i_H,a]
    if RKTLBM:
       for a in range(naR):
           f_r[0:my_i,0:mz_i,aR[a],0       ] = edgeL[i_F_r,a]
           f_b[0:my_i,0:mz_i,aR[a],0       ] = edgeL[i_F_b,a]
           g  [0:my_i,0:mz_i,aR[a],0       ] = edgeL[i_G,a]
           if rheology == "archean2":
              h[0:my_i,0:mz_i,aR[a],0      ] = edgeL[i_H,a]
    else: 
       for a in range(naR):
           f[0:my_i,0:mz_i,aR[a],0       ] = edgeL[i_F,a]
           g[0:my_i,0:mz_i,aR[a],0       ] = edgeL[i_G,a]
           if rheology == "archean2":
              h[0:my_i,0:mz_i,aR[a],0       ] = edgeL[i_H,a]

#-------------------------------------------------------------------------------
# Kernel to get the Z edges of a domain for overlap = False
#-------------------------------------------------------------------------------

def get_edgesZ():

    global i_F, i_G, i_H, i_F_r, i_F_b, i_HD, i_HB
    global i_rho_r, i_rho_b, i_D_ed, i_D_eb
    global    aL,    aR,    aU,    aD,    aB,    aF
    global   naL,   naR,   naU,   naD,   naB,   naF
    global edgeL    , edgeR    , edgeU    , edgeD    , edgeB    , edgeF
    global edge_rho_L    , edge_rho_R    , edge_rho_U    , edge_rho_D    , edge_rho_B    , edge_rho_F
    global   fgL    ,   fgR    ,   fgU    ,   fgD    ,   fgB    ,   fgF
    global   rhoL    ,  rhoR    ,  rhoU    ,  rhoD    ,  rhoB    ,  rhoF
    global R, Rx, Ry, Rz, rr, rl, rf, rb, rd, ru
    global reqSendR, reqSendL, recRecvR, reqReqvL
    global reqSendD, reqSendU, recRecvD, reqReqvU
    global reqSendF, reqSendB, recRecvF, reqReqvB

    if RKTLBM:
       for a in range(naU):
           fgD[i_F_r,a] = f_r[0:my_i,1     ,aU[a],0:mx_i]
           fgD[i_F_b,a] = f_b[0:my_i,1     ,aU[a],0:mx_i]
           fgD[i_G,a]   = g  [0:my_i,1     ,aU[a],0:mx_i]
           if rheology == "archean2":
              fgD[i_H,a] = h [0:my_i,1     ,aU[a],0:mx_i]
    else:
       for a in range(naU):
           fgD[i_F,a] = f[0:my_i,1     ,aU[a],0:mx_i]
           fgD[i_G,a] = g[0:my_i,1     ,aU[a],0:mx_i]
           if rheology == "archean2":
              fgD[i_H,a] = h[0:my_i,1     ,aU[a],0:mx_i]
    if RKTLBM:
       for a in range(naD):
           fgU[i_F_r,a] = f_r[0:my_i,mz[Rz],aD[a],0:mx_i]
           fgU[i_F_b,a] = f_b[0:my_i,mz[Rz],aD[a],0:mx_i]
           fgU[i_G,a]   = g  [0:my_i,mz[Rz],aD[a],0:mx_i]
           if rheology == "archean2":
              fgU[i_H,a] = h [0:my_i,mz[Rz],aD[a],0:mx_i]
    else:
       for a in range(naD):
           fgU[i_F,a] = f[0:my_i,mz[Rz],aD[a],0:mx_i]
           fgU[i_G,a] = g[0:my_i,mz[Rz],aD[a],0:mx_i]
           if rheology == "archean2":
              fgU[i_H,a] = h[0:my_i,mz[Rz],aD[a],0:mx_i]

    MPI.Prequest.Startall([reqSendD, reqRecvD, reqSendU, reqRecvU])
    MPI.Request.Waitall  ([reqSendD, reqRecvD, reqSendU, reqRecvU])
 
    if RKTLBM:
       for a in range(naU):
           f_r[0:my_i,mz[Rz]+1,aU[a],0:mx_i] = edgeD[i_F_r,a]
           f_b[0:my_i,mz[Rz]+1,aU[a],0:mx_i] = edgeD[i_F_b,a]
           g  [0:my_i,mz[Rz]+1,aU[a],0:mx_i] = edgeD[i_G,a]           
           if rheology == "archean2":
              h[0:my_i,mz[Rz]+1,aU[a],0:mx_i] = edgeD[i_H,a]
    else:
       for a in range(naU):
           f[0:my_i,mz[Rz]+1,aU[a],0:mx_i] = edgeD[i_F,a]
           g[0:my_i,mz[Rz]+1,aU[a],0:mx_i] = edgeD[i_G,a]
           if rheology == "archean2":
              h[0:my_i,mz[Rz]+1,aU[a],0:mx_i] = edgeD[i_H,a]     

    if RKTLBM:
       for a in range(naD):
           f_r[0:my_i,0       ,aD[a],0:mx_i] = edgeU[i_F_r,a]
           f_b[0:my_i,0       ,aD[a],0:mx_i] = edgeU[i_F_b,a]
           g  [0:my_i,0       ,aD[a],0:mx_i] = edgeU[i_G,a]
           if rheology == "archean2":
              h[0:my_i,0      ,aD[a],0:mx_i] = edgeU[i_H,a]
    else:
       for a in range(naD):
           f[0:my_i,0       ,aD[a],0:mx_i] = edgeU[i_F,a]
           g[0:my_i,0       ,aD[a],0:mx_i] = edgeU[i_G,a]
           if rheology == "archean2":
              h[0:my_i,0       ,aD[a],0:mx_i] = edgeU[i_H,a]
        
#-------------------------------------------------------------------------------
# Kernel to get the Y edges of a domain for overlap = False
#-------------------------------------------------------------------------------

def get_edgesY():

    global i_F, i_G, i_H, i_F_r, i_F_b, i_HD, i_HB
    global i_rho_r, i_rho_b, i_D_ed, i_D_eb
    global    aL,    aR,    aU,    aD,    aB,    aF
    global   naL,   naR,   naU,   naD,   naB,   naF
    global edgeL    , edgeR    , edgeU    , edgeD    , edgeB    , edgeF
    global edge_rho_L    , edge_rho_R    , edge_rho_U    , edge_rho_D    , edge_rho_B    , edge_rho_F
    global   fgL    ,   fgR    ,   fgU    ,   fgD    ,   fgB    ,   fgF
    global   rhoL    ,  rhoR    ,  rhoU    ,  rhoD    ,  rhoB    ,  rhoF
    global R, Rx, Ry, Rz, rr, rl, rf, rb, rd, ru
    global reqSendR, reqSendL, reqRecvR, reqRecvL
    global reqSendD, reqSendU, reqRecvD, reqRecvU
    global reqSendF, reqSendB, reqRecvF, reqRecvB

    if RKTLBM: 
       for a in range(naB):
           fgF[i_F_r,a] = f_r[1     ,0:mz_i,aB[a],0:mx_i]
           fgF[i_F_b,a] = f_b[1     ,0:mz_i,aB[a],0:mx_i]
           fgF[i_G,a]   = g  [1     ,0:mz_i,aB[a],0:mx_i]
           if rheology == "archean2":
              fgF[i_H,a] = h [1     ,0:mz_i,aB[a],0:mx_i]
    else:
       for a in range(naB):
           fgF[i_F,a]    = f[1     ,0:mz_i,aB[a],0:mx_i]
           fgF[i_G,a]    = g[1     ,0:mz_i,aB[a],0:mx_i]
           if rheology  == "archean2":
              fgF[i_H,a] = h[1  ,0:mz_i,aB[a],0:mx_i]
        
    if RKTLBM:
       for a in range(naF):
           fgB[i_F_r,a]  = f_r[my[Ry],0:mz_i,aF[a],0:mx_i]
           fgB[i_F_b,a]  = f_b[my[Ry],0:mz_i,aF[a],0:mx_i]
           fgB[i_G,a]    = g  [my[Ry],0:mz_i,aF[a],0:mx_i]
           if rheology == "archean2":
              fgB[i_H,a] = h  [my[Ry],0:mz_i,aF[a],0:mx_i]
    else: 
       for a in range(naF):
           fgB[i_F,a] = f[my[Ry],0:mz_i,aF[a],0:mx_i]
           fgB[i_G,a] = g[my[Ry],0:mz_i,aF[a],0:mx_i]
           if rheology == "archean2":
              fgB[i_H,a] = h[my[Ry],0:mz_i,aF[a],0:mx_i]

    MPI.Prequest.Startall([reqSendF, reqRecvF, reqSendB, reqRecvB])
    MPI.Request.Waitall  ([reqSendF, reqRecvF, reqSendB, reqRecvB])

    if RKTLBM:
       for a in range(naB):
           f_r[my[Ry]+1,0:mz_i,aB[a],0:mx_i] = edgeF[i_F_r,a]
           F_b[my[Ry]+1,0:mz_i,aB[a],0:mx_i] = edgeF[i_F_b,a]
           g  [my[Ry]+1,0:mz_i,aB[a],0:mx_i] = edgeF[i_G,a]
           if rheology == "archean2":
              h[my[Ry]+1,0:mz_i,aB[a],0:mx_i] = edgeF[i_H,a]
    else:
       for a in range(naB):
           f[my[Ry]+1,0:mz_i,aB[a],0:mx_i] = edgeF[i_F,a]
           g[my[Ry]+1,0:mz_i,aB[a],0:mx_i] = edgeF[i_G,a]
           if rheology == "archean2":
              h[my[Ry]+1,0:mz_i,aB[a],0:mx_i] = edgeF[i_H,a]
    if RKTLBM:
       for a in range(naF):
           f_r[0       ,0:mz_i,aF[a],0:mx_i] = edgeB[i_F_r,a]
           f_b[0       ,0:mz_i,aF[a],0:mx_i] = edgeB[i_F_b,a]
           g  [0       ,0:mz_i,aF[a],0:mx_i] = edgeB[i_G,a]
           if rheology == "archean2":
              h[0      ,0:mz_i,aF[a],0:mx_i] = edgeB[i_H,a]
    else:
       for a in range(naF):
           f[0       ,0:mz_i,aF[a],0:mx_i] = edgeB[i_F,a]
           g[0       ,0:mz_i,aF[a],0:mx_i] = edgeB[i_G,a]
           if rheology == "archean2":
              h[0       ,0:mz_i,aF[a],0:mx_i] = edgeB[i_H,a]
        

#-------------------------------------------------------------------------------
# Start send/recv requests to get edges for overlap = True
#-------------------------------------------------------------------------------

def get_edges_hd():
    # Standalone X-halo exchange for the advected-depletion populations hd.
    # Covers the x-direction only, so melt runs with hd must use
    # domain_decomposition = 1D. Column 1 goes to the left neighbour's right
    # halo and column mx to the right neighbour's left halo, matching the
    # packed f/g/h exchange (fgR -> rl, fgL -> rr). rl/rr are MPI.PROC_NULL
    # at non-periodic domain edges, which makes the Sendrecv a no-op there.
    mxx   = mxi   # per-rank physical width: hd.shape[3]-2 is wrong on the last rank (mx < nxi)
    sendL = np.ascontiguousarray(hd[:,:,:,1])
    sendR = np.ascontiguousarray(hd[:,:,:,mxx])
    recvR = np.empty_like(sendL)
    recvL = np.empty_like(sendR)
    comm.Sendrecv(sendL, dest=rl, sendtag=21, recvbuf=recvR, source=rr, recvtag=21)
    comm.Sendrecv(sendR, dest=rr, sendtag=22, recvbuf=recvL, source=rl, recvtag=22)
    hd[:,:,:,mxx+1] = recvR
    hd[:,:,:,0]     = recvL
    sendL = np.ascontiguousarray(hb[:,:,:,1])
    sendR = np.ascontiguousarray(hb[:,:,:,mxx])
    comm.Sendrecv(sendL, dest=rl, sendtag=23, recvbuf=recvR, source=rr, recvtag=23)
    comm.Sendrecv(sendR, dest=rr, sendtag=24, recvbuf=recvL, source=rl, recvtag=24)
    hb[:,:,:,mxx+1] = recvR
    hb[:,:,:,0]     = recvL

def get_edges_start(): # Start sending faces without face corners
                       # to overlap with interior calculations

    global i_F, i_G, i_H, i_F_r, i_F_b, i_HD, i_HB
    global R, Rx, Ry, Rz
    global xedge, yedge, zedge, nx_i, ny_i, nz_i, mxi, myi, mzi, mx_i, my_i, mz_i

#   Corners of cuboid

    global    edgeRD,    edgeRU,    edgeLD,    edgeLU
    global      fgRD,      fgRU,      fgLD,      fgLU
    global      naRD,      naRU,      naLD,      naLU
    global rrd,  aRD, rru,  aRU, rld,  aLD, rlu,  aLU
    global reqSendRD, reqSendRU, reqSendLD, reqSendLU
    global reqRecvRD, reqRecvRU, reqRecvLD, reqRecvLU

    global    edgeFD,    edgeFU,    edgeBD,    edgeBU
    global      fgFD,      fgFU,      fgBD,      fgBU
    global      naFD,      naFU,      naBD,      naBU
    global rfd,  aFD, rfu,  aFU, rbd,  aBD, rbu,  aBU
    global reqSendFD, reqSendFU, reqSendBD, reqSendBU
    global reqRecvFD, reqRecvFU, reqRecvBD, reqRecvBU

    global    edgeRF,    edgeRB,    edgeLF,    edgeLB
    global      fgRF,      fgRB,      fgLF,      fgLB
    global      naRF,      naRB,      naLF,      naLB
    global rrf,  aRF, rrb,  aRB, rlf,  aLF, rlb,  aLB
    global reqSendRF, reqSendRB, reqSendLF, reqSendLB
    global reqRecvRF, reqRecvRB, reqRecvLF, reqRecvLB

#   Save corner fg's ready to send

#   XZ corners
    if RKTLBM:
       for a in range(naLU):
           y1 = yedge ; y2 = myi+yedge
           fgRD[i_F_r] = f_r[y1:y2,1    ,aLU,1    ]
           fgLD[i_F_r] = f_r[y1:y2,1    ,aRU,mxi  ]
           fgRU[i_F_r] = f_r[y1:y2,mzi  ,aLD,1    ]
           fgLU[i_F_r] = f_r[y1:y2,mzi  ,aRD,mxi  ]
           fgRD[i_F_b] = f_b[y1:y2,1    ,aLU,1    ]
           fgLD[i_F_b] = f_b[y1:y2,1    ,aRU,mxi  ]
           fgRU[i_F_b] = f_b[y1:y2,mzi  ,aLD,1    ]
           fgLU[i_F_b] = f_b[y1:y2,mzi  ,aRD,mxi  ]
           fgRD[i_G] = g    [y1:y2,1    ,aLU,1    ]
           fgLD[i_G] = g    [y1:y2,1    ,aRU,mxi  ]
           fgRU[i_G] = g    [y1:y2,mzi  ,aLD,1    ]
           fgLU[i_G] = g    [y1:y2,mzi  ,aRD,mxi  ]
           if rheology == "archean2":
              fgRD[i_H] = h [y1:y2,1    ,aLU,1    ]
              fgLD[i_H] = h [y1:y2,1    ,aRU,mxi  ]
              fgRU[i_H] = h [y1:y2,mzi  ,aLD,1    ]
              fgLU[i_H] = h [y1:y2,mzi  ,aRD,mxi  ]
    else:
       for a in range(naLU):
           y1 = yedge ; y2 = myi+yedge
           fgRD[i_F] = f[y1:y2,1    ,aLU,1    ]
           fgRD[i_G] = g[y1:y2,1    ,aLU,1    ]
           fgLD[i_F] = f[y1:y2,1    ,aRU,mxi  ]
           fgLD[i_G] = g[y1:y2,1    ,aRU,mxi  ]
           fgRU[i_F] = f[y1:y2,mzi  ,aLD,1    ]
           fgRU[i_G] = g[y1:y2,mzi  ,aLD,1    ]
           fgLU[i_F] = f[y1:y2,mzi  ,aRD,mxi  ]
           fgLU[i_G] = g[y1:y2,mzi  ,aRD,mxi  ]
           if rheology == "archean2":
              fgRD[i_H] = h[y1:y2,1    ,aLU,1    ]
              fgLD[i_H] = h[y1:y2,1    ,aRU,mxi  ]
              fgRU[i_H] = h[y1:y2,mzi  ,aLD,1    ]
              fgLU[i_H] = h[y1:y2,mzi  ,aRD,mxi  ]

#   YZ corners
    if RKTLBM:
       for a in range(naBU):
        x1 = xedge ; x2 = mxi+xedge
        fgFD[i_F_r] = f_r[1    ,1    ,aBU,x1:x2]
        fgBD[i_F_r] = f_r[myi  ,1    ,aFU,x1:x2]
        fgFU[i_F_r] = f_r[1    ,mzi  ,aBD,x1:x2]
        fgBU[i_F_r] = f_r[myi  ,mzi  ,aFD,x1:x2]
        fgFD[i_F_b] = f_b[1    ,1    ,aBU,x1:x2]
        fgBD[i_F_b] = f_b[myi  ,1    ,aFU,x1:x2]
        fgFU[i_F_b] = f_b[1    ,mzi  ,aBD,x1:x2]
        fgBU[i_F_b] = f_b[myi  ,mzi  ,aFD,x1:x2]
        fgFD[i_G] = g    [1    ,1    ,aBU,x1:x2]
        fgBD[i_G] = g    [myi  ,1    ,aFU,x1:x2]
        fgFU[i_G] = g    [1    ,mzi  ,aBD,x1:x2]
        fgBU[i_G] = g    [myi  ,mzi  ,aFD,x1:x2]
        if rheology == "archean2":
           fgFD[i_H] = h [1    ,1    ,aBU,x1:x2]
           fgBD[i_H] = h [myi  ,1    ,aFU,x1:x2]
           fgFU[i_H] = h [1    ,mzi  ,aBD,x1:x2]
           fgBU[i_H] = h [myi  ,mzi  ,aFD,x1:x2]


    else:
       for a in range(naBU):
           x1 = xedge ; x2 = mxi+xedge
           fgFD[i_F] = f[1    ,1    ,aBU,x1:x2]
           fgFD[i_G] = g[1    ,1    ,aBU,x1:x2]
           fgBD[i_F] = f[myi  ,1    ,aFU,x1:x2]
           fgBD[i_G] = g[myi  ,1    ,aFU,x1:x2]
           fgFU[i_F] = f[1    ,mzi  ,aBD,x1:x2]
           fgFU[i_G] = g[1    ,mzi  ,aBD,x1:x2]
           fgBU[i_F] = f[myi  ,mzi  ,aFD,x1:x2]
           fgBU[i_G] = g[myi  ,mzi  ,aFD,x1:x2]
           if rheology == "archean2":
              fgFD[i_H] = h[1    ,1    ,aBU,x1:x2]
              fgBD[i_H] = h[myi  ,1    ,aFU,x1:x2]
              fgFU[i_H] = h[1    ,mzi  ,aBD,x1:x2]
              fgBU[i_H] = h[myi  ,mzi  ,aFD,x1:x2]

#   XY corners
    if RKTLBM:
       for a in range(naLB):
           z1 = zedge ; z2 = mzi+zedge
           fgRF[i_F_r] = f_r[1    ,z1:z2,aLB,1    ]
           fgLF[i_F_r] = f_r[1    ,z1:z2,aRB,mxi  ]
           fgRB[i_F_r] = f_r[myi  ,z1:z2,aLF,1    ]
           fgLB[i_F_r] = f_r[myi  ,z1:z2,aRF,mxi  ]
           fgRF[i_F_b] = f_b[1    ,z1:z2,aLB,1    ]
           fgLF[i_F_b] = f_b[1    ,z1:z2,aRB,mxi  ]
           fgRB[i_F_b] = f_b[myi  ,z1:z2,aLF,1    ]
           fgLB[i_F_b] = f_b[myi  ,z1:z2,aRF,mxi  ]
           fgRF[i_G] = g    [1    ,z1:z2,aLB,1    ]
           fgLF[i_G] = g    [1    ,z1:z2,aRB,mxi  ]
           fgRB[i_G] = g    [myi  ,z1:z2,aLF,1    ]
           fgLB[i_G] = g    [myi  ,z1:z2,aRF,mxi  ]
           if rheology == "archean2":
              fgRF[i_H] = h [1    ,z1:z2,aLB,1    ]
              fgLF[i_H] = h [1    ,z1:z2,aRB,mxi  ]
              fgRB[i_H] = h [myi  ,z1:z2,aLF,1    ]
              fgLB[i_H] = h [myi  ,z1:z2,aRF,mxi  ]           

    else:
       for a in range(naLB):
           z1 = zedge ; z2 = mzi+zedge
           fgRF[i_F] = f[1    ,z1:z2,aLB,1    ]
           fgRF[i_G] = g[1    ,z1:z2,aLB,1    ]
           fgLF[i_F] = f[1    ,z1:z2,aRB,mxi  ]
           fgLF[i_G] = g[1    ,z1:z2,aRB,mxi  ]
           fgRB[i_F] = f[myi  ,z1:z2,aLF,1    ]
           fgRB[i_G] = g[myi  ,z1:z2,aLF,1    ]
           fgLB[i_F] = f[myi  ,z1:z2,aRF,mxi  ]
           fgLB[i_G] = g[myi  ,z1:z2,aRF,mxi  ]
           if rheology == "archean2":
              fgRF[i_H] = h[1    ,z1:z2,aLB,1    ]
              fgLF[i_H] = h[1    ,z1:z2,aRB,mxi  ]
              fgRB[i_H] = h[myi  ,z1:z2,aLF,1    ]
              fgLB[i_H] = h[myi  ,z1:z2,aRF,mxi  ]
        
    if nrx > 1 and rank < nr:
        get_edgesX_start()   # Start sending X edges
    if nrz > 1 and rank < nr:
        get_edgesZ_start()   # Start sending Z edges
    if nry > 1 and rank < nr:
        get_edgesY_start()   # Start sending Y edges

#-------------------------------------------------------------------------------
# End send/recv requests to get edges for overlap = True
#-------------------------------------------------------------------------------

def get_edges_end():  # End sending cuboid edge faces after interior calculations
                      # and get cuboid corners of faces before  edge calculations

    global f_last, f, g_last, g, h_last, h, f_r_last, f_r, f_b_last, f_b

    if nrx > 1 and rank < nr:
        get_edgesX_end()    # Get X edges
    if nrz > 1 and rank < nr:
        get_edgesZ_end()    # Get Z edges
    if nry > 1 and rank < nr:
        get_edgesY_end()    # Get Y edges

    if rank < nr:
        get_corners()       # Get the corners

#-------------------------------------------------------------------------------
# Kernel to start send/recv requests of X edges of a domain for overlap = True
#-------------------------------------------------------------------------------

def get_edgesX_start():

    global i_F, i_G, i_H, i_F_r, i_F_b, i_HD, i_HB
    global    aL,    aR,    aU,    aD,    aB,    aF
    global   naL,   naR,   naU,   naD,   naB,   naF
    global edgeL    , edgeR    , edgeU    , edgeD    , edgeB    , edgeF
    global   fgL    ,   fgR    ,   fgU    ,   fgD    ,   fgB    ,   fgF
    global R, Rx, Ry, Rz, rr, rl, rf, rb, rd, ru
    global reqSendR, reqSendL, reqRecvR, reqRecvL
    global reqSendD, reqSendU, reqRecvD, reqRecvU
    global reqSendF, reqSendB, reqRecvF, reqRecvB
    global xedge, yedge, zedge, nx_i, ny_i, nz_i, mxi, myi, mzi, mx_i, my_i, mz_i

    if RKTLBM:
       for a in range(naL):
           fgR[i_F_r,a] = f_r[yedge:my_i-yedge,zedge:mz_i-zedge,aL[a],1     ]
           fgR[i_F_b,a] = f_b[yedge:my_i-yedge,zedge:mz_i-zedge,aL[a],1     ]
           fgR[i_G,a]   = g  [yedge:my_i-yedge,zedge:mz_i-zedge,aL[a],1     ]
           if rheology == "archean2":
              fgR[i_H,a] = h[yedge:my_i-yedge,zedge:mz_i-zedge,aL[a],1     ]
              
    else:
       for a in range(naL):
           fgR[i_F,a] = f[yedge:my_i-yedge,zedge:mz_i-zedge,aL[a],1     ]
           fgR[i_G,a] = g[yedge:my_i-yedge,zedge:mz_i-zedge,aL[a],1     ]
           if rheology == "archean2":
              fgR[i_H,a] = h[yedge:my_i-yedge,zedge:mz_i-zedge,aL[a],1     ]

    if RKTLBM:
       for a in range(naR):
           fgL[i_F_r,a] = f_r[yedge:my_i-yedge,zedge:mz_i-zedge,aR[a],mx[Rx]]
           fgL[i_F_b,a] = f_b[yedge:my_i-yedge,zedge:mz_i-zedge,aR[a],mx[Rx]]
           fgL[i_G,a]   = g  [yedge:my_i-yedge,zedge:mz_i-zedge,aR[a],mx[Rx]]
           if rheology == "archean2":
              fgL[i_H,a] = h [yedge:my_i-yedge,zedge:mz_i-zedge,aR[a],mx[Rx]]
 
    else:
       for a in range(naR):
           fgL[i_F,a] = f[yedge:my_i-yedge,zedge:mz_i-zedge,aR[a],mx[Rx]]
           fgL[i_G,a] = g[yedge:my_i-yedge,zedge:mz_i-zedge,aR[a],mx[Rx]]
           if rheology == "archean2":
              fgL[i_H,a] = h[yedge:my_i-yedge,zedge:mz_i-zedge,aR[a],mx[Rx]]
 

    MPI.Prequest.Startall([reqSendR, reqRecvR, reqSendL, reqRecvL])

#-------------------------------------------------------------------------------
# Kernel to end getting of X edges of a domain for overlap = True
#-------------------------------------------------------------------------------

def get_edgesX_end():

    global i_F, i_G, i_H, i_F_r, i_F_b, i_HD, i_HB
    global    aL,    aR,    aU,    aD,    aB,    aF
    global   naL,   naR,   naU,   naD,   naB,   naF
    global edgeL    , edgeR    , edgeU    , edgeD    , edgeB    , edgeF
    global   fgL    ,   fgR    ,   fgU    ,   fgD    ,   fgB    ,   fgF
    global R, Rx, Ry, Rz, rr, rl, rf, rb, rd, ru
    global reqSendR, reqSendL, reqRecvR, reqRecvL
    global reqSendD, reqSendU, reqRecvD, reqRecvU
    global reqSendF, reqSendB, reqRecvF, reqRecvB
    global xedge, yedge, zedge, nx_i, ny_i, nz_i, mxi, myi, mzi, mx_i, my_i, mz_i

    MPI.Request.Waitall  ([reqSendR, reqRecvR, reqSendL, reqRecvL])

    if RKTLBM:
       for a in range(naL):
           f_r_last[yedge:my_i-yedge,zedge:mz_i-zedge,aL[a],mx[Rx]+1] = edgeR[i_F_r,a]
           f_b_last[yedge:my_i-yedge,zedge:mz_i-zedge,aL[a],mx[Rx]+1] = edgeR[i_F_b,a]
           g_last  [yedge:my_i-yedge,zedge:mz_i-zedge,aL[a],mx[Rx]+1] = edgeR[i_G,a]
           if rheology == "archean2":
              h_last[yedge:my_i-yedge,zedge:mz_i-zedge,aL[a],mx[Rx]+1] = edgeR[i_H,a]

    else:
       for a in range(naL):
           f_last[yedge:my_i-yedge,zedge:mz_i-zedge,aL[a],mx[Rx]+1] = edgeR[i_F,a]
           g_last[yedge:my_i-yedge,zedge:mz_i-zedge,aL[a],mx[Rx]+1] = edgeR[i_G,a]
           if rheology == "archean2":
              h_last[yedge:my_i-yedge,zedge:mz_i-zedge,aL[a],mx[Rx]+1] = edgeR[i_H,a]
    if RKTLBM:
       for a in range(naR):
           f_r_last[yedge:my_i-yedge,zedge:mz_i-zedge,aR[a],0       ] = edgeL[i_F_r,a]
           f_b_last[yedge:my_i-yedge,zedge:mz_i-zedge,aR[a],0       ] = edgeL[i_F_b,a]
           g_last  [yedge:my_i-yedge,zedge:mz_i-zedge,aR[a],0       ] = edgeL[i_G,a]
           if rheology == "archean2":
              h_last[yedge:my_i-yedge,zedge:mz_i-zedge,aR[a],0       ] = edgeL[i_H,a]
    else:
       for a in range(naR):
           f_last[yedge:my_i-yedge,zedge:mz_i-zedge,aR[a],0       ] = edgeL[i_F,a]
           g_last[yedge:my_i-yedge,zedge:mz_i-zedge,aR[a],0       ] = edgeL[i_G,a]
           if rheology == "archean2":
              h_last[yedge:my_i-yedge,zedge:mz_i-zedge,aR[a],0       ] = edgeL[i_H,a]

#-------------------------------------------------------------------------------
# Kernel to start send/recv requests of Z edges of a domain for overlap = True
#-------------------------------------------------------------------------------

def get_edgesZ_start():

    global i_F, i_G, i_H, i_F_r, i_F_b, i_HD, i_HB
    global    aL,    aR,    aU,    aD,    aB,    aF
    global   naL,   naR,   naU,   naD,   naB,   naF
    global edgeL    , edgeR    , edgeU    , edgeD    , edgeB    , edgeF
    global   fgL    ,   fgR    ,   fgU    ,   fgD    ,   fgB    ,   fgF
    global R, Rx, Ry, Rz, rr, rl, rf, rb, rd, ru
    global reqSendR, reqSendL, recRecvR, reqReqvL
    global reqSendD, reqSendU, recRecvD, reqReqvU
    global reqSendF, reqSendB, recRecvF, reqReqvB
    global xedge, yedge, zedge, nx_i, ny_i, nz_i, mxi, myi, mzi, mx_i, my_i, mz_i

    if RKTLBM:
       for a in range(naU):
           fgD[i_F_r,a] = f_r[yedge:my_i-yedge,1     ,aU[a],xedge:mx_i-xedge]
           fgD[i_F_b,a] = f_b[yedge:my_i-yedge,1     ,aU[a],xedge:mx_i-xedge]
           fgD[i_G,a]   = g  [yedge:my_i-yedge,1     ,aU[a],xedge:mx_i-xedge]
           if rheology == "archean2":
              fgD[i_H,a] = h[yedge:my_i-yedge,1     ,aU[a],xedge:mx_i-xedge]
    else:
       for a in range(naU):
           fgD[i_F,a] = f[yedge:my_i-yedge,1     ,aU[a],xedge:mx_i-xedge]
           fgD[i_G,a] = g[yedge:my_i-yedge,1     ,aU[a],xedge:mx_i-xedge]
           if rheology == "archean2":
              fgD[i_H,a] = h[yedge:my_i-yedge,1     ,aU[a],xedge:mx_i-xedge]
    if RKTLBM:
       for a in range(naD):
           fgU[i_F_r,a] = f_r[yedge:my_i-yedge,mz[Rz],aD[a],xedge:mx_i-xedge]
           fgU[i_F_b,a] = f_b[yedge:my_i-yedge,mz[Rz],aD[a],xedge:mx_i-xedge]
           fgU[i_G,a]   = g  [yedge:my_i-yedge,mz[Rz],aD[a],xedge:mx_i-xedge]
           if rheology == "archean2":
              fgU[i_H,a] = h[yedge:my_i-yedge,mz[Rz],aD[a],xedge:mx_i-xedge]
    else:
       for a in range(naD):
           fgU[i_F,a] = f[yedge:my_i-yedge,mz[Rz],aD[a],xedge:mx_i-xedge]
           fgU[i_G,a] = g[yedge:my_i-yedge,mz[Rz],aD[a],xedge:mx_i-xedge]
           if rheology == "archean2":
              fgU[i_H,a] = h[yedge:my_i-yedge,mz[Rz],aD[a],xedge:mx_i-xedge]


    MPI.Prequest.Startall([reqSendD, reqRecvD, reqSendU, reqRecvU])

#-------------------------------------------------------------------------------
# Kernel to end getting of Z edges of a domain for overlap = True
#-------------------------------------------------------------------------------

def get_edgesZ_end():

    global i_F, i_G, i_H, i_F_r, i_F_b, i_HD, i_HB
    global    aL,    aR,    aU,    aD,    aB,    aF
    global   naL,   naR,   naU,   naD,   naB,   naF
    global edgeL    , edgeR    , edgeU    , edgeD    , edgeB    , edgeF
    global   fgL    ,   fgR    ,   fgU    ,   fgD    ,   fgB    ,   fgF
    global R, Rx, Ry, Rz, rr, rl, rf, rb, rd, ru
    global reqSendR, reqSendL, recRecvR, reqReqvL
    global reqSendD, reqSendU, recRecvD, reqReqvU
    global reqSendF, reqSendB, recRecvF, reqReqvB
    global xedge, yedge, zedge, nx_i, ny_i, nz_i, mxi, myi, mzi, mx_i, my_i, mz_i

    MPI.Request.Waitall  ([reqSendD, reqRecvD, reqSendU, reqRecvU])

    if RKTLBM:
       for a in range(naU):
           f_r_last[yedge:my_i-yedge,mz[Rz]+1,aU[a],xedge:mx_i-xedge] = edgeD[i_F_r,a]
           f_b_last[yedge:my_i-yedge,mz[Rz]+1,aU[a],xedge:mx_i-xedge] = edgeD[i_F_b,a]
           g_last  [yedge:my_i-yedge,mz[Rz]+1,aU[a],xedge:mx_i-xedge] = edgeD[i_G,a]
           if rheology == "archean2":
               h_last[yedge:my_i-yedge,mz[Rz]+1,aU[a],xedge:mx_i-xedge] = edgeD[i_H,a]
    else:
       for a in range(naU):
           f_last[yedge:my_i-yedge,mz[Rz]+1,aU[a],xedge:mx_i-xedge] = edgeD[i_F,a]
           g_last[yedge:my_i-yedge,mz[Rz]+1,aU[a],xedge:mx_i-xedge] = edgeD[i_G,a]
           if rheology == "archean2":
              h_last[yedge:my_i-yedge,mz[Rz]+1,aU[a],xedge:mx_i-xedge] = edgeD[i_H,a]

    if RKTLBM:
       for a in range(naD):
           f_r_last[yedge:my_i-yedge,0       ,aD[a],xedge:mx_i-xedge] = edgeU[i_F_r,a]
           f_b_last[yedge:my_i-yedge,0       ,aD[a],xedge:mx_i-xedge] = edgeU[i_F_b,a]
           g_last  [yedge:my_i-yedge,0       ,aD[a],xedge:mx_i-xedge] = edgeU[i_G,a]
           if rheology == "archean2":
              h_last[yedge:my_i-yedge,0       ,aD[a],xedge:mx_i-xedge] = edgeU[i_H,a]
    else:
       for a in range(naD):
           f_last[yedge:my_i-yedge,0       ,aD[a],xedge:mx_i-xedge] = edgeU[i_F,a]
           g_last[yedge:my_i-yedge,0       ,aD[a],xedge:mx_i-xedge] = edgeU[i_G,a]
           if rheology == "archean2":
              h_last[yedge:my_i-yedge,0       ,aD[a],xedge:mx_i-xedge] = edgeU[i_H,a]


#-------------------------------------------------------------------------------
# Kernel to start send/recv requests of Y edges of a domain for overlap = True
#-------------------------------------------------------------------------------

def get_edgesY_start():

    global i_F, i_G, i_H, i_F_r, i_F_b, i_HD, i_HB
    global    aL,    aR,    aU,    aD,    aB,    aF
    global   naL,   naR,   naU,   naD,   naB,   naF
    global edgeL    , edgeR    , edgeU    , edgeD    , edgeB    , edgeF
    global   fgL    ,   fgR    ,   fgU    ,   fgD    ,   fgB    ,   fgF
    global R, Rx, Ry, Rz, rr, rl, rf, rb, rd, ru
    global reqSendR, reqSendL, reqRecvR, reqRecvL
    global reqSendD, reqSendU, reqRecvD, reqRecvU
    global reqSendF, reqSendB, reqRecvF, reqRecvB
    global xedge, yedge, zedge, nx_i, ny_i, nz_i, mxi, myi, mzi, mx_i, my_i, mz_i
    

    if RKTLBM:
       for a in range(naB):
           fgF[i_F_r,a] = f_r[1     ,zedge:mz_i-zedge,aB[a],xedge:mx_i-xedge]
           fgF[i_F_b,a] = f_b[1     ,zedge:mz_i-zedge,aB[a],xedge:mx_i-xedge]
           fgF[i_G,a]   = g  [1     ,zedge:mz_i-zedge,aB[a],xedge:mx_i-xedge]
           if rheology == "archean2":
              fgF[i_H,a] = h [1     ,zedge:mz_i-zedge,aB[a],xedge:mx_i-xedge]
   
    else:
       for a in range(naB):
           fgF[i_F,a] = f[1     ,zedge:mz_i-zedge,aB[a],xedge:mx_i-xedge]
           fgF[i_G,a] = g[1     ,zedge:mz_i-zedge,aB[a],xedge:mx_i-xedge]
           if rheology == "archean2":
              fgF[i_H,a] = h[1     ,zedge:mz_i-zedge,aB[a],xedge:mx_i-xedge]
    
    if RKTLBM:
       for a in range(naF):
           fgB[i_F_r,a] = f_r[my[Ry],zedge:mz_i-zedge,aF[a],xedge:mx_i-xedge]
           fgB[i_F_b,a] = f_b[my[Ry],zedge:mz_i-zedge,aF[a],xedge:mx_i-xedge]
           fgB[i_G,a]   = g  [my[Ry],zedge:mz_i-zedge,aF[a],xedge:mx_i-xedge]
           if rheology == "archean2":
              fgB[i_H,a] = h[my[Ry],zedge:mz_i-zedge,aF[a],xedge:mx_i-xedge] 

    else:
       for a in range(naF):
           fgB[i_F,a] = f[my[Ry],zedge:mz_i-zedge,aF[a],xedge:mx_i-xedge]
           fgB[i_G,a] = g[my[Ry],zedge:mz_i-zedge,aF[a],xedge:mx_i-xedge]
           if rheology == "archean2":
              fgB[i_H,a] = h[my[Ry],zedge:mz_i-zedge,aF[a],xedge:mx_i-xedge]
 
    MPI.Prequest.Startall([reqSendF, reqRecvF, reqSendB, reqRecvB])

#-------------------------------------------------------------------------------
# Kernel to end getting of Y edges of a domain for overlap = True
#-------------------------------------------------------------------------------

def get_edgesY_end():

    global i_F, i_G, i_H, i_F_r, i_F_b, i_HD, i_HB
    global    aL,    aR,    aU,    aD,    aB,    aF
    global   naL,   naR,   naU,   naD,   naB,   naF
    global edgeL    , edgeR    , edgeU    , edgeD    , edgeB    , edgeF
    global   fgL    ,   fgR    ,   fgU    ,   fgD    ,   fgB    ,   fgF
    global R, Rx, Ry, Rz, rr, rl, rf, rb, rd, ru
    global reqSendR, reqSendL, reqRecvR, reqRecvL
    global reqSendD, reqSendU, reqRecvD, reqRecvU
    global reqSendF, reqSendB, reqRecvF, reqRecvB
    global xedge, yedge, zedge, nx_i, ny_i, nz_i, mxi, myi, mzi, mx_i, my_i, mz_i

    MPI.Request.Waitall  ([reqSendF, reqRecvF, reqSendB, reqRecvB])
    
    if RKTLBM:
       for a in range(naB):
           f_r_last[my[Ry]+1,zedge:mz_i-zedge,aB[a],xedge:mx_i-xedge] = edgeF[i_F_r,a]
           f_b_last[my[Ry]+1,zedge:mz_i-zedge,aB[a],xedge:mx_i-xedge] = edgeF[i_F_b,a]
           g_last  [my[Ry]+1,zedge:mz_i-zedge,aB[a],xedge:mx_i-xedge] = edgeF[i_G,a]
           if rheology == "archean2":
              h_last[my[Ry]+1,zedge:mz_i-zedge,aB[a],xedge:mx_i-xedge] = edgeF[i_H,a]
    else:
       for a in range(naB):
           f_last[my[Ry]+1,zedge:mz_i-zedge,aB[a],xedge:mx_i-xedge] = edgeF[i_F,a]
           g_last[my[Ry]+1,zedge:mz_i-zedge,aB[a],xedge:mx_i-xedge] = edgeF[i_G,a]
           if rheology == "archean2":
              h_last[my[Ry]+1,zedge:mz_i-zedge,aB[a],xedge:mx_i-xedge] = edgeF[i_H,a]

    if RKTLBM:
       for a in range(naF):
           f_r_last[0       ,zedge:mz_i-zedge,aF[a],xedge:mx_i-xedge] = edgeB[i_F_r,a]
           f_b_last[0       ,zedge:mz_i-zedge,aF[a],xedge:mx_i-xedge] = edgeB[i_F_b,a]
           g_lastvv[0       ,zedge:mz_i-zedge,aF[a],xedge:mx_i-xedge] = edgeB[i_G,a]
           if rheology == "archean2":
              h_last[0       ,zedge:mz_i-zedge,aF[a],xedge:mx_i-xedge] = edgeB[i_H,a]
    else:    
       for a in range(naF):
           f_last[0       ,zedge:mz_i-zedge,aF[a],xedge:mx_i-xedge] = edgeB[i_F,a]
           g_last[0       ,zedge:mz_i-zedge,aF[a],xedge:mx_i-xedge] = edgeB[i_G,a]
           if rheology == "archean2":
              h_last[0       ,zedge:mz_i-zedge,aF[a],xedge:mx_i-xedge] = edgeB[i_H,a]


#-------------------------------------------------------------------------------
# Kernel to get the corners of cuboid faces of domains (edges of cuboid)
# for overlap = True
#-------------------------------------------------------------------------------

def get_corners():

    global i_F, i_G, i_H, i_F_r, i_F_b, i_HD, i_HB
    global R, Rx, Ry, Rz
    global xedge, yedge, zedge, nx_i, ny_i, nz_i, mxi, myi, mzi, mx_i, my_i, mz_i

#   Corners of cuboid

    global    edgeRD,    edgeRU,    edgeLD,    edgeLU
    global      fgRD,      fgRU,      fgLD,      fgLU
    global      naRD,      naRU,      naLD,      naLU
    global rrd,  aRD, rru,  aRU, rld,  aLD, rlu,  aLU
    global reqSendRD, reqSendRU, reqSendLD, reqSendLU
    global reqRecvRD, reqRecvRU, reqRecvLD, reqRecvLU

    global    edgeFD,    edgeFU,    edgeBD,    edgeBU
    global      fgFD,      fgFU,      fgBD,      fgBU
    global      naFD,      naFU,      naBD,      naBU
    global rfd,  aFD, rfu,  aFU, rbd,  aBD, rbu,  aBU
    global reqSendFD, reqSendFU, reqSendBD, reqSendBU
    global reqRecvFD, reqRecvFU, reqRecvBD, reqRecvBU

    global    edgeRF,    edgeRB,    edgeLF,    edgeLB
    global      fgRF,      fgRB,      fgLF,      fgLB
    global      naRF,      naRB,      naLF,      naLB
    global rrf,  aRF, rrb,  aRB, rlf,  aLF, rlb,  aLB
    global reqSendRF, reqSendRB, reqSendLF, reqSendLB
    global reqRecvRF, reqRecvRB, reqRecvLF, reqRecvLB

    if naRD == 1:
        MPI.Prequest.Startall([reqSendLU, reqSendRU, reqSendLD, reqSendRD,
                               reqRecvLU, reqRecvRU, reqRecvLD, reqRecvRD])
        MPI.Prequest.Waitall ([reqSendLU, reqSendRU, reqSendLD, reqSendRD,
                               reqRecvLU, reqRecvRU, reqRecvLD, reqRecvRD])
    if naFD == 1:
        MPI.Prequest.Startall([reqSendBU, reqSendFU, reqSendBD, reqSendFD,
                               reqRecvBU, reqRecvFU, reqRecvBD, reqRecvFD])
        MPI.Prequest.Waitall ([reqSendBU, reqSendFU, reqSendBD, reqSendFD,
                               reqRecvBU, reqRecvFU, reqRecvBD, reqRecvFD])
    if naRF == 1:
        MPI.Prequest.Startall([reqSendLB, reqSendRB, reqSendLF, reqSendRF,
                               reqRecvLB, reqRecvRB, reqRecvLF, reqRecvRF])
        MPI.Prequest.Waitall ([reqSendLB, reqSendRB, reqSendLF, reqSendRF,
                               reqRecvLB, reqRecvRB, reqRecvLF, reqRecvRF])
#   XZ corners
    if RKTLBM:
       for a in range(naLU):
        y1 = yedge ; y2 = myi+yedge
        f_r_last[y1:y2,0    ,aRD,0    ] = edgeLU[i_F_r]
        f_r_last[y1:y2,mzi+1,aRU,0    ] = edgeLD[i_F_r]
        f_r_last[y1:y2,0    ,aLD,mxi+1] = edgeRU[i_F_r]
        f_r_last[y1:y2,mzi+1,aLU,mxi+1] = edgeRD[i_F_r]
        f_b_last[y1:y2,0    ,aRD,0    ] = edgeLU[i_F_b]
        f_b_last[y1:y2,mzi+1,aRU,0    ] = edgeLD[i_F_b]
        f_b_last[y1:y2,0    ,aLD,mxi+1] = edgeRU[i_F_b]
        f_b_last[y1:y2,mzi+1,aLU,mxi+1] = edgeRD[i_F_b]
        g_last  [y1:y2,0    ,aRD,0    ] = edgeLU[i_G]
        g_last  [y1:y2,mzi+1,aRU,0    ] = edgeLD[i_G]
        g_last  [y1:y2,0    ,aLD,mxi+1] = edgeRU[i_G]
        g_last  [y1:y2,mzi+1,aLU,mxi+1] = edgeRD[i_G]
        if rheology == "archean2":
           h_last[y1:y2,0    ,aRD,0    ] = edgeLU[i_H]
           h_last[y1:y2,mzi+1,aRU,0    ] = edgeLD[i_H]
           h_last[y1:y2,0    ,aLD,mxi+1] = edgeRU[i_H]
           h_last[y1:y2,mzi+1,aLU,mxi+1] = edgeRD[i_H]   
    else: 
       for a in range(naLU):
           y1 = yedge ; y2 = myi+yedge
           f_last[y1:y2,0    ,aRD,0    ] = edgeLU[i_F]
           g_last[y1:y2,0    ,aRD,0    ] = edgeLU[i_G]
           f_last[y1:y2,mzi+1,aRU,0    ] = edgeLD[i_F]
           g_last[y1:y2,mzi+1,aRU,0    ] = edgeLD[i_G]
           f_last[y1:y2,0    ,aLD,mxi+1] = edgeRU[i_F]
           g_last[y1:y2,0    ,aLD,mxi+1] = edgeRU[i_G]
           f_last[y1:y2,mzi+1,aLU,mxi+1] = edgeRD[i_F]
           g_last[y1:y2,mzi+1,aLU,mxi+1] = edgeRD[i_G]
           if rheology == "archean2":
              h_last[y1:y2,0    ,aRD,0    ] = edgeLU[i_H]
              h_last[y1:y2,mzi+1,aRU,0    ] = edgeLD[i_H]
              h_last[y1:y2,0    ,aLD,mxi+1] = edgeRU[i_H]
              h_last[y1:y2,mzi+1,aLU,mxi+1] = edgeRD[i_H]     

#   YZ corners
    if RKTLBM:
       for a in range(naBU):
           x1 = xedge ; x2 = mxi+xedge
           f_r_last[0    ,0    ,aFD,x1:x2] = edgeBU[i_F_r]
           f_r_last[0    ,mzi+1,aFU,x1:x2] = edgeBD[i_F_r]
           f_r_last[myi+1,0    ,aBD,x1:x2] = edgeFU[i_F_r]
           f_r_last[myi+1,mzi+1,aBU,x1:x2] = edgeFD[i_F_r]
           f_b_last[0    ,0    ,aFD,x1:x2] = edgeBU[i_F_b]
           f_b_last[0    ,mzi+1,aFU,x1:x2] = edgeBD[i_F_b]
           f_b_last[myi+1,0    ,aBD,x1:x2] = edgeFU[i_F_b]
           f_b_last[myi+1,mzi+1,aBU,x1:x2] = edgeFD[i_F_b]
           g_last  [0    ,0    ,aFD,x1:x2] = edgeBU[i_G]
           g_last  [0    ,mzi+1,aFU,x1:x2] = edgeBD[i_G]
           g_last  [myi+1,0    ,aBD,x1:x2] = edgeFU[i_G]
           g_last  [myi+1,mzi+1,aBU,x1:x2] = edgeFD[i_G]
           if rheology == "archean2":
              h_last[0    ,0    ,aFD,x1:x2] = edgeBU[i_H]
              h_last[0    ,mzi+1,aFU,x1:x2] = edgeBD[i_H]
              h_last[myi+1,0    ,aBD,x1:x2] = edgeFU[i_H]
              h_last[myi+1,mzi+1,aBU,x1:x2] = edgeFD[i_H]
    else:
       for a in range(naBU):
           x1 = xedge ; x2 = mxi+xedge
           f_last[0    ,0    ,aFD,x1:x2] = edgeBU[i_F]
           g_last[0    ,0    ,aFD,x1:x2] = edgeBU[i_G]
           f_last[0    ,mzi+1,aFU,x1:x2] = edgeBD[i_F]
           g_last[0    ,mzi+1,aFU,x1:x2] = edgeBD[i_G]
           f_last[myi+1,0    ,aBD,x1:x2] = edgeFU[i_F]
           g_last[myi+1,0    ,aBD,x1:x2] = edgeFU[i_G]
           f_last[myi+1,mzi+1,aBU,x1:x2] = edgeFD[i_F]
           g_last[myi+1,mzi+1,aBU,x1:x2] = edgeFD[i_G]
           if rheology == "archean2":
              h_last[0    ,0    ,aFD,x1:x2] = edgeBU[i_H]
              h_last[0    ,mzi+1,aFU,x1:x2] = edgeBD[i_H]        
              h_last[myi+1,0    ,aBD,x1:x2] = edgeFU[i_H]
              h_last[myi+1,mzi+1,aBU,x1:x2] = edgeFD[i_H]

#   XY corners
    if RKTLBM:
       for a in range(naLB):
           z1 = zedge ; z2 = mzi+zedge
           f_r_last[0    ,z1:z2,aRF,0    ] = edgeLB[i_F_r]
           f_r_last[myi+1,z1:z2,aRB,0    ] = edgeLF[i_F_r]
           f_r_last[0    ,z1:z2,aLF,mxi+1] = edgeRB[i_F_r]
           f_r_last[myi+1,z1:z2,aLB,mxi+1] = edgeRF[i_F_r]
           f_b_last[0    ,z1:z2,aRF,0    ] = edgeLB[i_F_b]
           f_b_last[myi+1,z1:z2,aRB,0    ] = edgeLF[i_F_b]
           f_b_last[0    ,z1:z2,aLF,mxi+1] = edgeRB[i_F_b]
           f_b_last[myi+1,z1:z2,aLB,mxi+1] = edgeRF[i_F_b]
           g_last  [0    ,z1:z2,aRF,0    ] = edgeLB[i_G]
           g_last  [myi+1,z1:z2,aRB,0    ] = edgeLF[i_G]
           g_last  [0    ,z1:z2,aLF,mxi+1] = edgeRB[i_G]
           g_last  [myi+1,z1:z2,aLB,mxi+1] = edgeRF[i_G]
           if rheology == "archean2":
              h_last[0    ,z1:z2,aRF,0    ] = edgeLB[i_H]
              h_last[myi+1,z1:z2,aRB,0    ] = edgeLF[i_H]
              h_last[0    ,z1:z2,aLF,mxi+1] = edgeRB[i_H]
              h_last[myi+1,z1:z2,aLB,mxi+1] = edgeRF[i_H]
    else:
       for a in range(naLB):
           z1 = zedge ; z2 = mzi+zedge
           f_last[0    ,z1:z2,aRF,0    ] = edgeLB[i_F]
           g_last[0    ,z1:z2,aRF,0    ] = edgeLB[i_G]
           f_last[myi+1,z1:z2,aRB,0    ] = edgeLF[i_F]
           g_last[myi+1,z1:z2,aRB,0    ] = edgeLF[i_G]
           f_last[0    ,z1:z2,aLF,mxi+1] = edgeRB[i_F]
           g_last[0    ,z1:z2,aLF,mxi+1] = edgeRB[i_G]
           f_last[myi+1,z1:z2,aLB,mxi+1] = edgeRF[i_F]
           g_last[myi+1,z1:z2,aLB,mxi+1] = edgeRF[i_G]
           if rheology == "archean2":
              h_last[0    ,z1:z2,aRF,0    ] = edgeLB[i_H]       
              h_last[myi+1,z1:z2,aRB,0    ] = edgeLF[i_H]
              h_last[0    ,z1:z2,aLF,mxi+1] = edgeRB[i_H]
              h_last[myi+1,z1:z2,aLB,mxi+1] = edgeRF[i_H]

#=======================================================================================
# Get the edges of the density fields onto the regions in each node
# This is needed for the RK LBM due to colour gradient calculations
#=======================================================================================

#-------------------------------------------------------------------------------
# get_edges_rho
#-------------------------------------------------------------------------------

def get_edges_rho():

    if nrx > 1 and rank < nr:
        get_edges_rho_X()   # Get rhoX edges
    if nrz > 1 and rank < nr:
        get_edges_rho_Z()   # Get rhoZ edges
    if nry > 1 and rank < nr:
        get_edges_rho_Y()   # Get rhoY edges

#-------------------------------------------------------------------------------
# Kernel to get the X edges of a domain
#-------------------------------------------------------------------------------

def get_edges_rho_X():

    global i_F_r, i_F_b, i_G
    global i_rho_r, i_rho_b, i_D_ed, i_D_eb
    global    aL,    aR,    aU,    aD,    aB,    aF
    global   naL,   naR,   naU,   naD,   naB,   naF
    global edgeL    , edgeR    , edgeU    , edgeD    , edgeB    , edgeF
    global edge_rho_L    , edge_rho_R    , edge_rho_U    , edge_rho_D    , edge_rho_B    , edge_rho_F
    global   fgL    ,   fgR    ,   fgU    ,   fgD    ,   fgB    ,   fgF
    global   rhoL    ,  rhoR    ,  rhoU    ,  rhoD    ,  rhoB    ,  rhoF
    global rho_r, rho_b, D_ed, D_eb, i_D_ed, i_D_eb
    global R, Rx, Ry, Rz, rr, rl, rf, rb, rd, ru
    global reqSendrhoR, reqSendrhoL, reqRecvrhoR, reqRecvrhoL
    global reqSendrhoD, reqSendrhoU, reqRecvrhoD, reqRecvrhoU
    global reqSendrhoF, reqSendrhoB, reqRecvrhoF, reqRecvrhoB

    rhoR[i_rho_r]   = rho_r[0:my_i,0:mz_i,1     ]
    rhoR[i_rho_b]   = rho_b[0:my_i,0:mz_i,1     ]
    rhoL[i_rho_r]   = rho_r[0:my_i,0:mz_i,mx[Rx]]
    rhoL[i_rho_b]   = rho_b[0:my_i,0:mz_i,mx[Rx]]
    if rheology == "archean2":
        rhoR[i_D_ed] = D_ed[0:my_i,0:mz_i,1     ]
        rhoR[i_D_eb] = D_eb[0:my_i,0:mz_i,1     ]
        rhoL[i_D_ed] = D_ed[0:my_i,0:mz_i,mx[Rx]]
        rhoL[i_D_eb] = D_eb[0:my_i,0:mz_i,mx[Rx]]

    MPI.Prequest.Startall([reqSendrhoR, reqRecvrhoR, reqSendrhoL, reqRecvrhoL])
    MPI.Request.Waitall  ([reqSendrhoR, reqRecvrhoR, reqSendrhoL, reqRecvrhoL])

    rho_r[0:my_i,0:mz_i,mx[Rx]+1] = edge_rho_R[i_rho_r,:,:]
    rho_b[0:my_i,0:mz_i,mx[Rx]+1] = edge_rho_R[i_rho_b,:,:]
    rho_r[0:my_i,0:mz_i,0       ] = edge_rho_L[i_rho_r,:,:]
    rho_b[0:my_i,0:mz_i,0       ] = edge_rho_L[i_rho_b,:,:]
    if rheology == "archean2":
        D_ed[0:my_i,0:mz_i,mx[Rx]+1] = edge_rho_R[i_D_ed,:,:]
        D_eb[0:my_i,0:mz_i,mx[Rx]+1] = edge_rho_R[i_D_eb,:,:]
        D_ed[0:my_i,0:mz_i,0       ] = edge_rho_L[i_D_ed,:,:]
        D_eb[0:my_i,0:mz_i,0       ] = edge_rho_L[i_D_eb,:,:]

#-------------------------------------------------------------------------------
# Kernel to get the Z edges of a domain
#-------------------------------------------------------------------------------

def get_edges_rho_Z():

    global i_F_r, i_F_b, i_G
    global i_rho_r, i_rho_b, i_D_ed, i_D_eb
    global    aL,    aR,    aU,    aD,    aB,    aF
    global   naL,   naR,   naU,   naD,   naB,   naF
    global edgeL    , edgeR    , edgeU    , edgeD    , edgeB    , edgeF
    global edge_rho_L    , edge_rho_R    , edge_rho_U    , edge_rho_D    , edge_rho_B    , edge_rho_F
    global   fgL    ,   fgR    ,   fgU    ,   fgD    ,   fgB    ,   fgF
    global   rhoL    ,  rhoR    ,  rhoU    ,  rhoD    ,  rhoB    ,  rhoF
    global R, Rx, Ry, Rz, rr, rl, rf, rb, rd, ru
    global reqSendrhoR, reqSendrhoL, reqRecvrhoR, reqRecvrhoL
    global reqSendrhoD, reqSendrhoU, reqRecvrhoD, reqRecvrhoU
    global reqSendrhoF, reqSendrhoB, reqRecvrhoF, reqRecvrhoB

    rhoD[i_rho_r]   = rho_r[0:my_i,1     ,0:mx_i]
    rhoD[i_rho_b]   = rho_b[0:my_i,1     ,0:mx_i]
    rhoU[i_rho_r]   = rho_r[0:my_i,mz[Rz],0:mx_i]
    rhoU[i_rho_b]   = rho_b[0:my_i,mz[Rz],0:mx_i]

    MPI.Prequest.Startall([reqSendrhoD, reqRecvrhoD, reqSendrhoU, reqRecvrhoU])
    MPI.Request.Waitall  ([reqSendrhoD, reqRecvrhoD, reqSendrhoU, reqRecvrhoU])

    rho_r[0:my_i,mz[Rz]+1,0:mx_i] = edge_rho_D[i_rho_r]
    rho_b[0:my_i,mz[Rz]+1,0:mx_i] = edge_rho_D[i_rho_b]
    rho_r[0:my_i,0       ,0:mx_i] = edge_rho_U[i_rho_r]
    rho_b[0:my_i,0       ,0:mx_i] = edge_rho_U[i_rho_b]

#-------------------------------------------------------------------------------
# Kernel to get the Y edges of a domain
#-------------------------------------------------------------------------------

def get_edges_rho_Y():

    global i_F_r, i_F_b, i_G
    global i_rho_r, i_rho_b, i_D_ed, i_D_eb
    global    aL,    aR,    aU,    aD,    aB,    aF
    global   naL,   naR,   naU,   naD,   naB,   naF
    global edgeL    , edgeR    , edgeU    , edgeD    , edgeB    , edgeF
    global edge_rho_L    , edge_rho_R    , edge_rho_U    , edge_rho_D    , edge_rho_B    , edge_rho_F
    global   fgL    ,   fgR    ,   fgU    ,   fgD    ,   fgB    ,   fgF
    global   rhoL    ,  rhoR    ,  rhoU    ,  rhoD    ,  rhoB    ,  rhoF
    global R, Rx, Ry, Rz, rr, rl, rf, rb, rd, ru
    global reqSendrhoR, reqSendrhoL, reqRecvrhoR, reqRecvrhoL
    global reqSendrhoD, reqSendrhoU, reqRecvrhoD, reqRecvrhoU
    global reqSendrhoF, reqSendrhoB, reqRecvrhoF, reqRecvrhoB

    rhoF[i_rho_r]   = rho_r[1     ,0:mz_i,0:mx_i]
    rhoF[i_rho_b]   = rho_b[1     ,0:mz_i,0:mx_i]
    rhoB[i_rho_r]   = rho_r[my[Ry],0:mz_i,0:mx_i]
    rhoB[i_rho_b]   = rho_b[my[Ry],0:mz_i,0:mx_i]

    MPI.Prequest.Startall([reqSendrhoF, reqRecvrhoF, reqSendrhoB, reqRecvrhoB])
    MPI.Request.Waitall  ([reqSendrhoF, reqRecvrhoF, reqSendrhoB, reqRecvrhoB])

    rho_r[my[Ry]+1,0:mz_i,0:mx_i] = edge_rho_F[i_rho_r]
    rho_b[my[Ry]+1,0:mz_i,0:mx_i] = edge_rho_F[i_rho_b]
    rho_r[0       ,0:mz_i,0:mx_i] = edge_rho_B[i_rho_r]
    rho_b[0       ,0:mz_i,0:mx_i] = edge_rho_B[i_rho_b]

#===============================================================================
# Functions for various rheologies or physics
#===============================================================================
#-------------------------------------------------------------------------------
# Heat from radiogenic decay
#-------------------------------------------------------------------------------

# USERS can change the heat from radiogenic decay in this function if required

def heat_from_radiogenic_decay(t):

#   Calculate the mass of one lattice site in kg and secs per computational time
#   step which allows physical constants to be converted into lattice units

    m_p = rho_p * dx_p**D          # kg (mass per lattice site)
    secs_per_time_step = dt_p*dt   # Secs/(computational time step = dt)

#   Heat capacity Cp from Fabio A. Capitanio et al., 2022, JGR: SE
#   and radiogenic heat generation rate H are used to calculate the
#   heat rate H_r in K/s. This can then be converted to H_r in lattice units
#   Are these numbers right? They come to only 480 degrees per Gy.

    H_r = H/Cp                      # J/(kg.s)/(J/(kg.K)) = K/s

#   Note that regarding of units, one can calculate Cp_lattice=Cp_phys*m_p (J/K) and
#   H_lattice = H_phys*m_p/secs_per_time_step (J/(lattice mass x lattice time step)
#   but these constants are not needed in lattice units as H_r converts to K/s which
#   can then be converted to K/(lattice time step) as below. H_r then be added
#   to g_alpha at each time step to allow for radiogenic heat generation.

    H_r *= secs_per_time_step       # K/(lattice time step)

    return np_float(H_r)

#=======================================================================================
# Main program
#
# Thermal LBM with option for temperature dependent viscosity and/or plastic rheology
#
# Authors:      Peter Mora and Gabriele Morra
#
# COPYRIGHT
#=======================================================================================

#--------------------------------------------------------------
# MPI setup
#--------------------------------------------------------------

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

#================================================================
# Get PROJECT and RUN names and make directory names

if len(argv) == 3:
    script, PROJECT, RUN = argv
else:

#   Self documentation option if not 2 arguments

    if rank == 0:
        HELP = "no"
        if len(argv) == 2:
            script, HELP = argv

        print(file=sys.stderr,flush=True)
        print("============================================================================",file=sys.stderr,flush=True)

        if HELP != "-h":
            print("ERROR: Incorrect usage of TLBM_fast.py",file=sys.stderr,flush=True)
            print("============================================================================",file=sys.stderr,flush=True)

        print("TLBM_fast.py documentation",file=sys.stderr,flush=True)
        print("============================================================================",file=sys.stderr,flush=True)
        print("Usage: mpiexec -np CORES python3 TLBM_fast.py PROJECT RUN  (do a run)",file=sys.stderr,flush=True)
        print("       mpiexec -np CORES python3 TLBM_fast.py -h          (print documentation)",file=sys.stderr,flush=True)
        print("----------------------------------------------------------------------------",file=sys.stderr,flush=True)
        print("Notes CORES   = number of cores",file=sys.stderr,flush=True)
        print("      PROJECT = project name",file=sys.stderr,flush=True)
        print("      RUN     = run name",file=sys.stderr,flush=True)
        print("----------------------------------------------------------------------------",file=sys.stderr,flush=True)
        print("Input  file: ./projects/PROJECT/RUN.in         (must exist)",file=sys.stderr,flush=True)
        print("Output file: ./projects/PROJECT/RUN.out        (will be created)",file=sys.stderr,flush=True)
        print("Data files:  ./projects/PROJECT/data/RUN*.npy  (will be created)",file=sys.stderr,flush=True)
        print("Plot files:  ./projects/PROJECT/plots/RUN*.png (will be created)",file=sys.stderr,flush=True)
        print("============================================================================",file=sys.stderr,flush=True)
    exit()

#==============================================================

# Set directory names

proj_dir  = "./projects/" + PROJECT + "/"
run_dir   = "./projects/" + PROJECT + "/" + RUN + "/"
movie_dir = run_dir + "movies/"
data_dir  = run_dir + "data/"
plot_dir  = run_dir + "plots/"

# Set the strings for the input and output file names

in_file    = proj_dir + RUN + ".in"
out_file   = proj_dir + RUN + ".out"
bmk_file   = proj_dir + RUN + ".bmk" # Benchmark data of w_TBL & V_rms
data_files = data_dir  + "*.npy"
plot_files = plot_dir + "*.png"

#==============================================================
# Check the directories exist and if necessary make them

#   Exit if the project input directory and input run file does not exist

if not os.path.exists(in_file):
    if rank == 0:
        print("ERROR: project input file for run does not exist",file=sys.stderr,flush=True)
        print("       You must create it to do a run",file=sys.stderr,flush=True)
    exit()

#   Make any required output directories

if rank == 0:
    if not os.path.exists(proj_dir):  os.mkdir(proj_dir)
    if not os.path.exists(run_dir):   os.mkdir(run_dir)
    if not os.path.exists(movie_dir): os.mkdir(movie_dir)
    if not os.path.exists(data_dir):  os.mkdir(data_dir)
    if not os.path.exists(plot_dir):  os.mkdir(plot_dir)

# Set stdout to the out file and print banner for run
# and open the benchmark file

if rank == 0:
    f_bmk = open(bmk_file,"w")       # Open  benchmark file
    f_bmk.writelines("# it t_My log10Ra dt T_bar Nu Nu_u Nu_d Vrms_nd v_rms_cm_py mach\n")
    f_bmk.close()
    sys.stdout = open(out_file, "w") # Set stdout to output file
    for iout in range(2):
        if   iout == 0: output_print = sys.stderr
        elif iout == 1: output_print = sys.stdout
        print("==================================================================",file=output_print,flush=True)
        print("TLBM_fast.py version %s (released %s), updated %s " %(release,release_date,version_date),file=output_print,flush=True)
        print("Starting %d core run" %(size),file=output_print,flush=True)
        now = datetime.now()
        date_string =  now.strftime("%B %d, %Y, %H:%M:%S")
        print("------------------------------------------------------------------",file=output_print,flush=True)
        print("Run date: %s" %(date_string),file=output_print,flush=True)

        print("Run = %s, project = %s" %(RUN,PROJECT),file=output_print,flush=True)
        print("==================================================================",file=output_print,flush=True)
        print("Input  file  = ",in_file,file=output_print,flush=True)
        print("Data   files = ",data_files,file=output_print,flush=True)
        print("Plot   files = ",plot_files,file=output_print,flush=True)
        print("Output file  = ",out_file,file=output_print,flush=True)
        print("==================================================================",file=output_print,flush=True)

#===============================================================================
# Set default specifications
#===============================================================================

#-------------------------------------------------------------------------------
# The specifications below are whole mantle values
#-------------------------------------------------------------------------------

dT_p       = 3400                 # Temperature difference         (K)
L_p        = 2900*1000            # Mantle thickness               (m)
rho_p      = 4500                 # Mantle density                 (kg/m^3)
rho_z_file = "unspecified"        # Density depth profile           (kg/m^3)
T_z_file   = "unspecified"        # Temperature depth profile       (C)
P_z_file   = "unspecified"        # Pressure depth profile          (Pa)
L_z_file       = "unspecified"    # Lithology depth distribution profile
basalt_depth_p = 0                # Initial basalt layer thickness at surface (m); 0 = no initial basalt
eta_p      = 2.8e22               # Dynamic viscosity              (Pa.s = kg/(m.s)
ka_p       = 0.5e-6               # Thermal diffusivity            (m^2/s)

#===============================================================================
# Specifications below are from the archean example
#===============================================================================

#-------------------------------------------------------------------------------
# Run title and description on plots
#-------------------------------------------------------------------------------

run_title            = "Archean"                               # Run title on plots
description          = "Archean_example_(author:_Peter_Mora)"  # Description of the run on plots (printed at lower right of image)
description_color    = "black"                                 # Color of the description
description_box      = False                                   # Plot the description in a white box
description_location = "frame"                                 # Description location (image or frame)

#-------------------------------------------------------------------------------
# Model size, timesteps and geometry
#-------------------------------------------------------------------------------

nt                   = "5000000"      # Number of timesteps, physical time will be timesteps x dt_p where dt_p = physical time per timestep
R_p                  = "Earth"        # Radius of the planet
L_p                  = "upper_mantle" # Thickness of the layer
nx                   = 512            # X resolution of model
ny                   = 1              # Y resolution of model
nz                   = 64             # Z resolution of model
geometry             = "annular"      # Geometry, cartesian or spherical (if ny=1, spherical means annular)
theta1               = 65             # Start theta of an annular model (degrees)
theta2               = 115            # End   theta of an annular model (degrees)
n_TBL                = 3              # Number of Thermal Boundary Layer (TBL) points to use to calculate the Nusselt Number

#-------------------------------------------------------------------------------
# Physical specifications
#-------------------------------------------------------------------------------

state                = "solid"
g_dx_gauge           = False         # Annular g~dx^3 buoyancy gauge (benchmark test_option S); single-phase path only       # State of the fluid for the simulation: "solid" or "magma"
power_correction     = 1.1           # Correction power for mode = "magma" to obtain velocities consitent with studies
Ra                   = 1.e5          # Rayleigh number: float or "unspecified"
Pr                   = 2000          # Prandtl number:  float
Peclet             = "unspecified"   # Peclet number
rho_p                = 3300          # Density                         (kg/m^3)
eta_p                = 5e21          # Dynamic viscosity               (Pa.s = kg/(m.s))
eta_min_p            = 1e18          # Viscosity lower clamp           (Pa.s)
eta_max_p            = 1e28          # Viscosity upper clamp           (Pa.s)
# Arrhenius viscosity profile (3 depth zones: upper mantle z<660km, lower mantle 660<=z<2740km, core region z>=2740km)
act_E_1              = 300e3         # Activation energy zone 1 (upper mantle)   (J/mol)
act_E_2              = 370e3         # Activation energy zone 2 (lower mantle)   (J/mol)
act_E_3              = 162e3         # Activation energy zone 3 (core region)    (J/mol)
act_V_1              = 5.00e-6       # Activation volume zone 1 (upper mantle)   (m^3/mol)
act_V_2              = 3.65e-6       # Activation volume zone 2 (lower mantle)   (m^3/mol)
act_V_3              = 1.40e-6       # Activation volume zone 3 (core region)    (m^3/mol)
Delta_eta_1          = 1.0           # Viscosity multiplier zone 1 (upper mantle)
Delta_eta_2          = 30.0          # Viscosity multiplier zone 2 (lower mantle)
Delta_eta_3          = 0.1           # Viscosity multiplier zone 3 (core region)
P_scale_1            = 1e15          # Pressure scale zone 1 (upper mantle)      (Pa)
P_scale_2            = 200e9         # Pressure scale zone 2 (lower mantle)      (Pa)
P_scale_3            = 1610e9        # Pressure scale zone 3 (core region)       (Pa)
ka_p               = "unspecified"   # Thermal diffusivity             (m^2/s)
g_p                  = 9.8           # Gravity                         (m/s^2)
beta_p               = 3e-5          # Thermal expansion coefficient   (1/K)
T0_p                 = 15            # Surface temperature             (C)
dT_p                 = 1800          # Temperature difference          (K)
dT_p_random          = 10            # Random temperature fluctuations (K) on thermal BC to initiate convection
T_p_random           = 0             # Random temperature fluctuations (K) throughout the model
initialize_ncells    = 0             # Initialize this many convection cells by sinusoidal fluctuation on thermal BC
cell_amplitude       = 1             # Amplitude of sinusoidal temperature fluctuations on thermal BC if initialize_ncells > 0
amp_sinusoid         = 0.0           # Amplitude of sinusoidal temperature fluctuations throughut the model
x_periods            = 1             # X periods of sinusoidal temperature fluctuations if amp_sinusoid > 0
y_periods            = 1             # Y periods of sinusoidal temperature fluctuations if amp_sinusoid > 0
H                    = 1.87e-11      # Radiogenic heat generation rate (W/kg = J/(kg.s))
Cp                   = 1200          # Heat capacity                   (J/(kg.K) (physical units))
temperature          = "reynolds"    # Temperature dependence of viscosity
gamma_T              = 10            # Power in reynolds temperature dependence
rheology             = "archean"     # Rheology
nu_formula           = "Tau"         # Option: "Tau", nu_plastic = (Tau_y/Tau)^N or "strain_rate", nu_plastic = nu_star + (Tau_y/(2 * strain_rate))^N
Tau_y0_p             = 5e7           # Yield stress                    (Pa = kg/(m.s^2))
Tau_yi_p             = 1e3           # Yield stress gradient           (Pa/m = kg/(m^2.s^2)), 1 MPa/km = 10^3 x Pa/m)
nu_star              = 0             #  
N                    = 4.0           # Exponent for plastic rheology equation: Tau = (Tau_y/Tau)^N
RKTLBM               = False         # Multiphase RKTLBM 
collision            = True          # 2nd collision term on/off
melt_time            = 120           # Beginning of melt extraction
melt_position        = "surface"     # Only if archean2. Where is melt added: 1) surface 2) insitu 3) off  
melt_extraction      = False         # Turm on and off melt
melt_test_mode       = False         # TEST: zero real M_ext, inject synthetic melt to verify MPI transport
rho_r_file          = "unspecified"  # RKTLBM first density file
rho_b_file          = "unspecified"  # RKTLBM second density file
L_file               = "unspecified"  # 2D lithology IC npy in the run dir (archean2): L field (ny,nz,nx), values 0..1
L_file               = "unspecified"  # 2D lithology IC npy in the run dir (archean2): L field (ny,nz,nx), values 0..1
rho1                 = 1.0           # RKTLBM red density 
rho2                 = 1.0           # RKTLBM blue density
rho1_p               = 2750          # RKTLBM red physical density
rho2_p               = 3300          # RKTLBM blue physical density
rho12                = 0.8           # RKTLBM red/blue physical density
theta_u              = 90            # 90 degrees on upper surface
theta_d              = 90            # 90 degrees on lower surface
A                    = 0
interface_thickness  = 0.75 
nu_r_factor          = 1
beta_g_factor        = 0             # No convection
g_comp_factor        = 0             # No compositional gravity     

#-------------------------------------------------------------------------------
# Archean2 parameters needed for TTG melt extract calculations (Requested parameters by TLBM_fast32 and TLBM_fast64)
#-------------------------------------------------------------------------------
P1_thresh            = 1.6
P2_thresh            = 10
M1_pd                = 0.00
M2_pd                = 0.05
M3_pd                = 0.50
M4_pd                = 0.60
M1_bs                = 0.02
M2_bs                = 0.04
TTG_DI_max           = 1e9         # Global cap (%) on cumulative TTG DI; basalt->TTG extraction
                                   # pauses while DI_bs >= this (2.4 = present-day crust/upper-mantle
                                   # DI, Gunawardana et al. 2020). 1e9 = no cap.
M3_bs                = 0.50
M4_bs                = 0.60
rho0_pd              = 3300  # Peridotite density
rho0_bs              = 3000  # Basalt density 
rho0_ttg             = 2750  # TTG density 
rho_melt_bs          = 2750  # Melt-bearing hydrated basalt
rho_melt_pd          = 3000  # Melt-bearing peridotite 

DI_bs                = []
bs_capped            = False
TTG_LP_frac          = 1.0         # Max share (fraction of DI_bs) the LP window may hold; 1 = off
TTG_MP_frac          = 1.0         # Max share for the MP window
TTG_HP_frac          = 1.0         # Max share for the HP window
fields_every_n_data  = 1           # write npy fields every Nth data event (DI/quota update runs every event)
lp_gatev             = 1.0         # proportional window gates (1 = fully open)
mp_gatev             = 1.0
hp_gatev             = 1.0
DI_bs_LP             = []
DI_bs_MP             = []
DI_bs_HP             = []
DI_pd                = []
Mex_LP_percent       = []
Mex_MP_percent       = []
Mex_HP_percent       = []

#-------------------------------------------------------------------------------
# Boundary conditions,  fixed, free_slip, or periodic (side boundaries only)
#-------------------------------------------------------------------------------

BC_sides             = "periodic"    # Side  BC
BC_upper             = "free_slip"   # Upper BC
BC_lower             = "fixed"       # Lower BC

#-------------------------------------------------------------------------------
# Physical timestep specifications
#-------------------------------------------------------------------------------

time_per_timestep    = "500_Years"   # Physical time per timestep (preferred, takes precedence over nu_l & g_l)
nu_l                 = "unspecified" # Lattice viscosity   (use nu_l to calculate g_l  & hence, dt_p, takes prescedence over g_l)
g_l                  = "unspecified" # Lattice gravity     (use g_l  to calculate nu_l & hence  dt_p)
factor_dt            = 1.0           # Factor to modify dt (only applied if time_per_timestep = unspecified)

#-------------------------------------------------------------------------------
# Mach specifications to maximize speed & control accuracy/stability
#-------------------------------------------------------------------------------

factor_max_vary_dt   = 30          # Maximum factor to vary dt
Inc_vary_dt          = "2000"      # Time step increment to vary dt to control mach
start_vary_dt        = "20000"     # Start time step when step dt can start to be varied
factor_increase_dt   = 1.05        # Factor increase in dt if mach is too low
factor_decrease_dt   = 0.75        # Factor decrease in dt if mach is too high
max_mach_limit       = 0.01        # Upper mach limit before dt will be decreased
min_mach_limit       = 0.005       # Lower mach limit before dt will be increased
mach_limit_to_exit   = 0.5         # Maximum mach number before a run will exit 
nu_l_min             = 0.00015     # Minimum allowable lattice viscosity

#-------------------------------------------------------------------------------
# Calculation and accuracy specifications
#-------------------------------------------------------------------------------

float_bits           = 32          # Precision in bits of calculations
edge_bits            = 32          # Precision in bits of edges of domains
data_bits            = 16          # Precision in bits of data files (.npy) of plots (.png)
shifted_distribution = True        # Increases LBM accuracy
pre_streaming_BC     = False       # Use pre-streaming wall BCs (matches RKTLBM strmac):
                                   #   True  = read wall-node f_last at the wall (correct
                                   #           halfway bounce-back; same as RKTLBM strmac).
                                   #   False = legacy: read post-streamed neighbour value.

#-------------------------------------------------------------------------------
# Output specifications
#-------------------------------------------------------------------------------

output_plots         = True        # Output plot files
output_data          = False       # Output data files
output_units         = "physical"  # Output units for prints & plots: physical or lattice
output_timing        = True        # Output performance & timing

#-------------------------------------------------------------------------------
# Output increment specifications in timesteps, physical increment is time_per_timestep x timesteps
#-------------------------------------------------------------------------------

Inc_plots            = "10000"     # Plot increment
Inc_data             = "10000"     # Data increment
Inc_timing           = "50000"     # Timing increment
Inc_mach             = "1000"      # Mach number print/check increment

#-------------------------------------------------------------------------------
# Output plot specifications
#-------------------------------------------------------------------------------

plane                = "xz"        # Plane to plot (xz, xy or yz)
coordinates          = "0"         # Coordinate (ie. if plane=xz, it is the y-coordinate of the plane plotted)

plot_size            = "small"     # Plot size: vsmall, small, medium, large, vlarge
plot_size_factor     = 1.0         # Plot size modifying factor to apply to plot_size
plot_colorbar        = True        # Plot the colorbar
colormap             = "jet"       # Colormap for plots: eg, jet, bwr, seismic, ...
latex                = False       # Use latex for titles and labels in plots

plot_rho             = False       # Plot the density     (kg/m^3)
plot_T               = True        # Plot the temperature (degrees C)
plot_u               = True        # Plot the magnitude of velocity  |u| (cm/year)
plot_nu_total        = False       # Plot total viscosity
plot_nu_kinematic    = False       # Plot kinematic viscosity
plot_nu_dynamic      = False       # Plot dynamic viscosity
plot_Ra              = False       # Plot per-cell Rayleigh number (log10)
plot_vectors         = False       # Plot velocity vectors
plot_Tau             = False       # Plot shear stress (only for code = python)
plot_L               = False       # Plot density energy distribution 
plot_rho_eff         = False       # Plot effective density (kg/m^3)
plot_rho_s           = False       # Plot solidus density   (kg/m^3)
plot_M_pd            = False       # Peridotite melt fraction
plot_M_bs            = False       # TTG  melt fraction 
plot_M_MP            = False       # MP-TTG melt 
plot_M_LP            = False       # LP-TTG melt
plot_M_HP            = False       # HP-TTG melt
plot_ce_LP           = False       # Cumulative extracted LP TTG
plot_ce_MP           = False       # Cumulative extracted MP TTG
plot_ce_HP           = False       # Cumulative extracted HP TTG
plot_ce_pd           = False       # Cumulative extracted peridotite
plot_de_pd           = False       # Advected depletion field D_ed (peridotite)
plot_de_bs           = False       # Advected depletion field D_eb (basalt)
plot_ce_bs           = False       # Cumulative extracted peridotite
plot_T_sol_pd        = False       # T_solidus peridotite
plot_T_sol_bs        = False       # T_solidus basalt 
plot_T_liq_pd        = False       # T_liquidus peridotite
plot_T_liq_bs        = False       # T_liquidus basalt
plot_P               = False       # Pressure
plot_DI              = False       # Plot Differentiation Index
plot_rho_r           = False       # Plot red density
plot_rho_b           = False       # Plot blue density
plot_f_r             = False       
plot_f_b             = False
plot_FF_mag          = False       # Plot color gradient
plot_FF_x            = False       # Plot color gradient
plot_FF_z            = False       # Plot color gradient
plot_FF_y            = False       # Plot color gradient

#-------------------------------------------------------------------------------
# Output plot scale specifications
#
# Plots show scaled fields based on the max and min values in the field,
# and two specifications, clip (linear scaling) and gpow (nonlinear scaling).
#
# clip: clips a field above clip times the maximum value in the field
# EXAMPLE: clip = 0.5 clips at half the maximum so values less than half
# of the maximum will have the full color range and values above half
# the maximum show as saturated colors.

# gpow: a nonlinear scaling of field f by a factor of (|f|/max(f))^gpow.
# EXAMPLE: gpow = 0.5 boosts small values by the square root of magnitude
# of the field allowing more detail to be seen of the smaller values
#-------------------------------------------------------------------------------

vector_overlay       = True        # Overlay the velocity vectors on |u| plot
vectors_color        = "w"         # Color of the overlay vectors on |u| plot
vector_length        = 1.0         # Length of vectors relative to the default
vector_width         = 1.0         # Width  of vectors relative to the default
vector_interval      = 6           # Interval in lattice units to plot vectors
annulus_interval     = 10          # Interval in annnular lattice plot

clip_rho             = 1.0         # clip for density
gpow_rho             = 1.0         # gpow for density
clip_T               = 1.0         # clip for temperature
gpow_T               = 1.0         # gpow for temperature
fix_scale_T          = True        # Fix the temperature scale or not
clip_u               = 1.0         # clip for the magnitude of velocity
gpow_u               = 1.0         # gpow for the magnitude of velocity
clip_nu_total        = 1.0         # clip for nu_total
gpow_nu_total        = 1.0         # gpow for nu_total
clip_Tau             = 1.0         # clip for Tau
gpow_Tau             = 1.0         # gpow for Tau

#-------------------------------------------------------------------------------
# Performance specifications
#-------------------------------------------------------------------------------

code                 = "c"           # Do calculations using c or python
overlap              = False         # Overlap getting of domain edges with calculations. This only
                                     # has potential benefit for large cluster runs with code = c
domain_decomposition = "3D"          # 1D, 2D or 3D. 3D will optimize the domain decomposition (DD)
nxi_min              = 256           # Minimum x domain scan in DD (min cache blocking size)
nxi_max              = 512           # Maximum x domain scan in DD (max cache blocking size)
nyi_min              = "unspecified" # Minimum y domain scan in DD (useful to specify for huge 3D runs as faster)
nyi_max              = "unspecified" # Maximum y domain scan in DD (useful to specify for huge 3D runs as faster)
nzi_min              = "unspecified" # Minimum z domain scan in DD (useful to specify for huge 3D runs as faster)
nzi_max              = "unspecified" # Maximum z domain scan in DD (useful to specify for huge 3D runs as faster)

#===============================================================================
# Get the run specifications from the input file and import necessary functions
#===============================================================================

import get_parameters
from get_parameters import get_parameters
from get_parameters import str2str
get_parameters(in_file,rank)
from get_parameters import *
coordinates = coordinates.split(",")

import utilities
from utilities import sec_ky
from utilities import sec_y
from utilities import sec_my
from utilities import u_cm_py
from utilities import timevalue_to_timesteps
from utilities import timetext
from utilities import timevalue
from utilities import timevalue_to_secs
from utilities import timesteps_to_physical_time_text
from utilities import time_to_physical_time_text
from utilities import get_plane
from utilities import get_lattice_viscosity
from utilities import convert_units
from utilities import get_u_convert
from utilities import get_annular_geometry
from utilities import get_R_and_L

#-------------------------------------------------------------------------------
# Check for any invalid specifications
#-------------------------------------------------------------------------------

fatal_input_error = False

error = False
nt           , timescale_nt           , error = timevalue("nt"           ,nt           ,error,rank)

Inc_vary_dt  , timescale_Inc_vary_dt  , error = timevalue("Inc_vary_dt"  ,Inc_vary_dt  ,error,rank)
start_vary_dt, timescale_start_vary_dt, error = timevalue("start_vary_dt",start_vary_dt,error,rank)
Inc_plots    , timescale_Inc_plots    , error = timevalue("Inc_plots"    ,Inc_plots    ,error,rank)
Inc_data     , timescale_Inc_data     , error = timevalue("Inc_data"     ,Inc_data     ,error,rank)
Inc_mach     , timescale_Inc_mach     , error = timevalue("Inc_mach"     ,Inc_mach     ,error,rank)
Inc_timing   , timescale_Inc_timing   , error = timevalue("Inc_timing"   ,Inc_timing   ,error,rank)

if error: fatal_input_error = True

for iout in range(2):
    if   iout == 0: output_print = sys.stderr
    elif iout == 1: output_print = sys.stdout
    if not (BC_sides == "free_slip" or BC_sides == "fixed" or BC_sides == "periodic"):
        fatal_input_error = True
        if rank == 0: print('ERROR: BC_sides must be either free_slip, fixed or periodic',file=output_print,flush=True)
    if not (BC_upper == "free_slip" or BC_upper == "fixed"):
        if rank == 0: print('ERROR: BC_upper must be either free_slip or fixed',file=output_print,flush=True)
        fatal_input_error = True
    if not (BC_lower == "free_slip" or BC_lower == "fixed"):
        if rank == 0: print('ERROR: BC_lower must be free_slip or fixed',file=output_print,flush=True)
        fatal_input_error = True

    if plane != "xz" and plane != "xy" and plane != "yz":
        if rank == 0: print('ERROR: plane must be xz, xy or yz',file=output_print,flush=True)
        fatal_input_error = True

    if not (float_bits == 32 or float_bits == 64):
        if rank == 0: print("ERROR: float_bits must be 32 or 64",file=output_print,flush=True)
        fatal_input_error = True

    if not (edge_bits == 16 or edge_bits == 32 or edge_bits == 64):
        if rank == 0: print("ERROR: edge_bits must be 16, 32 or 64",file=output_print,flush=True)
        fatal_input_error = True

    if not (data_bits == 16 or data_bits == 32 or data_bits == 64):
        if rank == 0: print("ERROR: data_bits must be 16 or 32",file=output_print,flush=True)
        fatal_input_error = True

    if state != "solid" and state != "magma":
        if rank == 0: print("ERROR: state must be specified",file=output_print,flush=True)
        fatal_input_error = True

    if Ra == "unspecified" and eta_p == "unspecified":
        if rank == 0: print("ERROR: If Ra = unspecified, then eta_p must be specified",file=output_print,flush=True)
        fatal_input_error = True

#    if Ra != "unspecified" and eta_p != "unspecified":
#        if rank == 0: print("ERROR: Only one of Ra and eta_p can be specified",file=output_print,flush=True)
#        fatal_input_error = True

    if not (temperature == "reynolds" or temperature == "arrhenius" or temperature == "unspecified"):
        if rank == 0: print('ERROR: temperature can only be reynolds or arrhenius or unspecified',file=output_print,flush=True)
        fatal_input_error = True

    if not (rheology == "linear" or rheology == "archean" or rheology == "archean2"):
        if rank == 0: print('ERROR: rheology can only be linear, archean or archean2',file=output_print,flush=True)
        fatal_input_error = True

    if (time_per_timestep == "unspecified" and g_l == "unspecified" and nu_l == "unspecified"):
        if rank == 0: print('ERROR: One of g_l, nu_l or time_per_timestep must be specified',file=output_print,flush=True)
        fatal_input_error = True

    if Ra == "unspecified" and ka_p == "unspecified":
        if rank == 0: print('ERROR: If Ra is unspecified, ka_p must be specified',file=output_print,flush=True)
        fatal_input_error = True

    if ka_p == "unspecified":
        if rank == 0: print('ERROR: ka_p must be specified',file=output_print,flush=True)
        fatal_input_error = True

    if Pr == "unspecified":
        if rank == 0: print('ERROR: Pr must be specified',file=output_print,flush=True)
        fatal_input_error = True

if fatal_input_error:
    for iout in range(2):
        if   iout == 0: output_print = sys.stderr
        elif iout == 1: output_print = sys.stdout
        if rank == 0:
            print("------------------------------------------------------------------",file=output_print,flush=True)
            print("The input file has errors which must be fixed to do a run",file=output_print,flush=True)
            now = datetime.now()
            date_string =  now.strftime("%B %d, %Y, %H:%M:%S")
            print("==================================================================",file=output_print,flush=True)
            print("The run has exited: %s" %(date_string),file=output_print,flush=True)
            print("==================================================================",file=output_print,flush=True)
            print(file=output_print,flush=True)
    exit()
else:
    for iout in range(2):
        if   iout == 0: output_print = sys.stderr
        elif iout == 1: output_print = sys.stdout
        if rank == 0:
            now = datetime.now()
            date_string =  now.strftime("%B %d, %Y, %H:%M:%S")
            print("Run is starting: %s" %(date_string),file=output_print,flush=True)
            print("==================================================================",file=output_print,flush=True)

grep_mpi    = subprocess.run(['ps', '-alF'], stdout=subprocess.PIPE, text=True)
grep_output = [line for line in grep_mpi.stdout.splitlines() if 'mpi' in line]
if rank == 0:
   print("------------------------------------------------------------------",file=sys.stderr,flush=True)
   for line in grep_output:
       print("MPI process is running: %s" %(line),file=sys.stderr,flush=True)
       print("MPI process is running: %s" %(line),file=output_print,flush=True)
   print("------------------------------------------------------------------",file=sys.stderr,flush=True)

#-------------------------------------------------------------------------------
# Set np data types for ints, floats, edges and data
#-------------------------------------------------------------------------------

np_int = np.int32

if   float_bits == 32: np_float = np.float32
elif float_bits == 64: np_float = np.float64

if   edge_bits  == 16: np_edge  = np.float16
elif edge_bits  == 32: np_edge  = np.float32
elif edge_bits  == 64: np_edge  = np.float64

if   data_bits  == 16: np_data = np.float16
elif data_bits  == 32: np_data = np.float32
elif data_bits  == 64: np_data = np.float64

#-------------------------------------------------------------------------------
# Set MPI types for send/recv of floats, edges, and data
#-------------------------------------------------------------------------------

if   float_bits == 32: MPI_FLOAT = MPI.FLOAT
elif float_bits == 64: MPI_FLOAT = MPI.DOUBLE

if   edge_bits  == 16: MPI_EDGE  = MPI.BYTE.Create_contiguous(2).Commit()
elif edge_bits  == 32: MPI_EDGE  = MPI.FLOAT
elif edge_bits  == 64: MPI_EDGE  = MPI.DOUBLE

if   data_bits  == 16: MPI_DATA  = MPI.BYTE.Create_contiguous(2).Commit()
elif data_bits  == 32: MPI_DATA  = MPI.FLOAT
elif data_bits  == 64: MPI_DATA  = MPI.DOUBLE

# End   of get input parameters
#===============================================================================
# Initializations of c function calls
#===============================================================================

# Array definitions

array_5d_float = np.ctypeslib.ndpointer(np_float,ndim=5,flags='aligned,contiguous')
array_4d_float = np.ctypeslib.ndpointer(np_float,ndim=4,flags='aligned,contiguous')
array_3d_float = np.ctypeslib.ndpointer(np_float,ndim=3,flags='aligned,contiguous')
array_2d_float = np.ctypeslib.ndpointer(np_float,ndim=2,flags='aligned,contiguous')
array_1d_float = np.ctypeslib.ndpointer(np_float,ndim=1,flags='aligned,contiguous')
array_5d_int   = np.ctypeslib.ndpointer(np_int  ,ndim=5,flags='aligned,contiguous')
array_4d_int   = np.ctypeslib.ndpointer(np_int  ,ndim=4,flags='aligned,contiguous')
array_3d_int   = np.ctypeslib.ndpointer(np_int  ,ndim=3,flags='aligned,contiguous')
array_2d_int   = np.ctypeslib.ndpointer(np_int  ,ndim=2,flags='aligned,contiguous')
array_1d_int   = np.ctypeslib.ndpointer(np_int,ndim=1,flags='aligned,contiguous')
array_4d_bool  = np.ctypeslib.ndpointer(bool,ndim=4,flags='aligned,contiguous')
array_3d_bool  = np.ctypeslib.ndpointer(bool,ndim=3,flags='aligned,contiguous')
array_4d_short = np.ctypeslib.ndpointer(np.int16,ndim=4,flags='aligned,contiguous')

if   float_bits == 32:
    C_FLOAT = c_float
elif float_bits == 64:
    C_FLOAT = c_double

# Load the shared object file of the C TLBM_fast function

TLBM_fast = CDLL('./TLBM_fast.so')

# The arguments below show the meaning of the arguments in the argtypes specifications
#
#  TLBM_fast(
#  na,D,nx_i,ny_i,nz_i,mx_i,my_i,mz_i,xedge,yedge,zedge, \
#  w,c,shifted_distribution,overlap,interior, \
#  rho,eta,u, \
#  f_last,f,g_last,g, \
#  x_xa,y_ya,z_za, \
#  nu_f, tau_g, tau_h, BETA_G, \
#  Temperature, Rheology, Nu_formula, \
#  nu_l_min, gamma_T, Tau_y, N, Delta_g, \
#  nu_f_last, tau_f_last, nu_total, \
#  nu_star, \
#  P1_thresh, P2_thresh, \
#  M1_pd, M2_pd, M3_pd, M4_pd, \
#  M1_bs, M2_bs, M3_bs, M4_bs, \
#  rho0_pd, rho_melt_pd, M_ext_pd, \
#  rho0_bs, rho_melt_bs, M_ext_bs, \
#  M_ext_bs_LP_TTG, M_ext_bs_MP_TTG, M_ext_bs_HP_TTG, \
#  mask_bs, mask_pd, \
#  g_l, \
#  T_nd, \
#  cumulative_extracted_pd, cumulative_extracted_bs, \
#  cumulative_extracted_LP, cumulative_extracted_MP, cumulative_extracted_HP,\
#  P0, rho_M, L_ed,\
#  h, h_last, \
#  DR_2, T0_p, dT_p, z_nd, \
#  rho_solid, rho_eff, \
#  boundary_l,xl,naL,nbL,bL,biL, \
#  boundary_r,xr,naR,nbR,bR,biR, \
#  boundary_b,yb,naB,nbB,bB,biB, \
#  boundary_f,yf,naF,nbF,bF,biF, \
#  boundary_u,zu,naU,nbU,bU,biU, \
#  boundary_d,zd,naD,nbD,bD,biD, \
#  eta_u,eta_d, \
#  rank
#  )

TLBM_fast_argtypes = [ \
    c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int, \
    array_1d_float, array_2d_int, c_bool, c_bool, c_bool, \
    array_3d_float, array_3d_float, array_4d_float, \
    array_4d_float, array_4d_float, array_4d_float, array_4d_float, \
    array_2d_int, array_2d_int, array_2d_int, \
    array_1d_float, array_1d_float, array_1d_float, array_1d_float, \
    c_char_p, c_char_p, c_char_p, \
    C_FLOAT, C_FLOAT, C_FLOAT, C_FLOAT, array_1d_float, array_1d_float, array_1d_float, \
    array_3d_float, array_3d_float, array_3d_float, \
    C_FLOAT, \
    C_FLOAT, C_FLOAT, \
    C_FLOAT, C_FLOAT, C_FLOAT, C_FLOAT, \
    C_FLOAT, C_FLOAT, C_FLOAT, C_FLOAT, \
    C_FLOAT, C_FLOAT, array_3d_float, \
    C_FLOAT, C_FLOAT, array_3d_float, \
    array_3d_float, array_3d_float, array_3d_float, \
    array_3d_bool, array_3d_bool, \
    C_FLOAT, \
    array_3d_float, \
    array_3d_float, array_3d_float, \
    array_3d_float, array_3d_float, array_3d_float,\
    array_3d_float, array_3d_float, array_3d_float,\
    array_4d_float, array_4d_float, \
    C_FLOAT, C_FLOAT, C_FLOAT, array_3d_float,   \
    array_3d_float, array_3d_float, \
    c_bool, c_int, c_int, c_int, array_1d_int, array_1d_int, \
    c_bool, c_int, c_int, c_int, array_1d_int, array_1d_int, \
    c_bool, c_int, c_int, c_int, array_1d_int, array_1d_int, \
    c_bool, c_int, c_int, c_int, array_1d_int, array_1d_int, \
    c_bool, c_int, c_int, c_int, array_1d_int, array_1d_int, \
    c_bool, c_int, c_int, c_int, array_1d_int, array_1d_int, \
    array_2d_float, array_2d_float, \
    c_bool, \
    c_int \
]

#  The arguments below show the meaning of the arguments in the RKTLBM_fast_strmac argtypes specifications
#
#  RKTLBM_strmac(
#  na,D,nx_i,ny_i,nz_i,xedge,yedge,zedge, \
#  mx_i,my_i,mz_i, \
#  g_last,g, \
#  x_xa,y_ya,z_za, \
#  rho_r, rho_b,
#  f_r_last, f_r, \
#  f_b_last, f_b, \
#  rheology,
#  h_last,h, \
#  boundary_l,xl,naL,nbL,bL,biL, \
#  boundary_r,xr,naR,nbR,bR,biR, \
#  boundary_b,yb,naB,nbB,bB,biB, \
#  boundary_f,yf,naF,nbF,bF,biF, \
#  boundary_u,zu,naU,nbU,bU,biU, \
#  boundary_d,zd,naD,nbD,bD,biD, \
#  rank
#  )

RKTLBM_fast_strmac_argtypes = [ \
    c_int,c_int,c_int,c_int,c_int,c_int,c_int,c_int,\
    c_int,c_int,c_int,\
    array_4d_float, array_4d_float,\
    array_2d_int,array_2d_int, array_2d_int,\
    array_3d_float, array_3d_float,\
    array_4d_float, array_4d_float,\
    array_4d_float, array_4d_float,\
    c_char_p,\
    array_4d_float, array_4d_float,\
    array_4d_float, array_4d_float,\
    array_4d_float, array_4d_float,\
    c_bool, c_int, c_int, c_int, array_1d_int, array_1d_int, \
    c_bool, c_int, c_int, c_int, array_1d_int, array_1d_int, \
    c_bool, c_int, c_int, c_int, array_1d_int, array_1d_int, \
    c_bool, c_int, c_int, c_int, array_1d_int, array_1d_int, \
    c_bool, c_int, c_int, c_int, array_1d_int, array_1d_int, \
    c_bool, c_int, c_int, c_int, array_1d_int, array_1d_int, \
    C_FLOAT, C_FLOAT, C_FLOAT, C_FLOAT, \
    c_bool, \
    array_3d_float, array_3d_float, \
    c_int\
]

# The arguments below show the meaning of the arguments in the RKTLBM_fast_rest argtypes specifications
#
#  RKTLBM_rest( \
#  na,D,nx_i,ny_i,nz_i,xedge,yedge,zedge, \
#  mx_i,my_i,mz_i, \
#  w, c, \
#  eta, \
#  u, \
#  g, \
#  x_xa,y_ya,z_za, \
#  beta_g, \
#  tau_f0_r, tau_f0_b, tau_g, \
#  rho_r, rho_b, \
#  f_r, \
#  f_b, \
#  C_r, C_b, \
#  rho1_p, rho2_p, rho3_p, beta_g_factor, g_comp_factor, \
#  FF, \
#  A, B, interface_thickness, g_comp, \
#  BC_r_wetting_u, BC_b_wetting_u, BC_r_wetting_d, BC_b_wetting_d, \
#  Rheology, Temperature, nu_formula\
#  tau_h, \
#  nu_f_r, nu_f_b, nu_star\
#  nu_l_min, gamma_T, Tau_y, N, Delta_g, \ 
#  nu_f_last, tau_f_last, \
#  T_nd, \
#  act_E, act_V, Delta_eta_array, eta_p_input, DR_2, T0_p, dT_p, z_nd, \
#  h, L_ed, \
#  P0, cumulative_extracted_pd, cumulative_extracted_bs, \
#  cumulative_extracted_LP, cumulative_extracted_MP, cumulative_extracted_HP, \
#  P1_thresh, P2_thresh, \
#  M1_pd, M2_pd, M4_pd, \
#  M1_bs, M2_bs, M4_bs, \
#  collision, mask_bs, mask_pd, \
#  M_ext_bs_LP_TTG, M_ext_bs_MP_TTG, M_ext_bs_HP_TTG,
#  rho0_pd, rho_melt_pd, M_ext_pd, \
#  rho0_bs, rho_melt_bs, \
#  rho_solid, rho_eff, \
#  g_l, rho_M, rho_reff, \
#  melt_extraction, Melt_position, \
#  boundary_l,xl,naL,nbL,bL,biL, \
#  boundary_r,xr,naR,nbR,bR,biR, \
#  boundary_b,yb,naB,nbB,bB,biB, \
#  boundary_f,yf,naF,nbF,bF,biF, \
#  boundary_u,zu,naU,nbU,bU,biU, \
#  boundary_d,zd,naD,nbD,bD,biD, \
#  eta_u,eta_d, \
#  rank)

RKTLBM_fast_rest_argtypes = [ \
    c_int,c_int,c_int,c_int,c_int,c_int,c_int,c_int,\
    c_int,c_int,c_int,\
    array_1d_float, array_2d_int,\
    array_3d_float,\
    array_4d_float, \
    array_4d_float, \
    array_2d_int,array_2d_int, array_2d_int,\
    C_FLOAT,\
    C_FLOAT, C_FLOAT, array_3d_float,\
    array_3d_float,array_3d_float,\
    array_4d_float,\
    array_4d_float,\
    array_1d_float,array_1d_float,\
    C_FLOAT,C_FLOAT, C_FLOAT, C_FLOAT, C_FLOAT, \
    array_4d_float,\
    C_FLOAT, array_1d_float, C_FLOAT, C_FLOAT,\
    C_FLOAT, C_FLOAT, C_FLOAT, C_FLOAT, \
    c_char_p,c_char_p, c_char_p,\
    array_1d_float, \
    array_1d_float, array_1d_float, C_FLOAT,\
    C_FLOAT, C_FLOAT, C_FLOAT, C_FLOAT, array_1d_float, array_1d_float, array_1d_float, \
    array_3d_float, array_3d_float, \
    array_3d_float, \
    array_1d_float, array_1d_float, array_1d_float, C_FLOAT, C_FLOAT, C_FLOAT, C_FLOAT, array_3d_float, \
    array_4d_float, array_3d_float, \
    array_4d_float, array_3d_float, \
    array_4d_float, array_3d_float, \
    array_3d_float, array_3d_float, array_3d_float, \
    array_3d_float, array_3d_float, array_3d_float, \
    C_FLOAT, C_FLOAT, \
    C_FLOAT, C_FLOAT, C_FLOAT,\
    C_FLOAT, C_FLOAT, C_FLOAT,\
    c_bool, array_3d_bool, array_3d_bool, \
    array_3d_float, array_3d_float, array_3d_float, \
    C_FLOAT, C_FLOAT, array_3d_float,\
    C_FLOAT, C_FLOAT, \
    array_3d_float, array_3d_float, \
    C_FLOAT, array_3d_float, array_3d_float, \
    c_bool, c_char_p, \
    c_bool, c_int, c_int, c_int, array_1d_int, array_1d_int, \
    c_bool, c_int, c_int, c_int, array_1d_int, array_1d_int, \
    c_bool, c_int, c_int, c_int, array_1d_int, array_1d_int, \
    c_bool, c_int, c_int, c_int, array_1d_int, array_1d_int, \
    c_bool, c_int, c_int, c_int, array_1d_int, array_1d_int, \
    c_bool, c_int, c_int, c_int, array_1d_int, array_1d_int, \
    array_2d_float, array_2d_float,\
    c_bool, \
    array_3d_float, array_3d_float, \
    c_int,\
    C_FLOAT, C_FLOAT, C_FLOAT\
]

if   float_bits == 32:
     if RKTLBM:
        TLBM_fast.RKTLBM_strmac32.argtypes = RKTLBM_fast_strmac_argtypes
        RKTLBM_strmac                      = TLBM_fast.RKTLBM_strmac32
        TLBM_fast.RKTLBM_rest32.argtypes   = RKTLBM_fast_rest_argtypes
        RKTLBM_rest                        = TLBM_fast.RKTLBM_rest32
     else:
        TLBM_fast.TLBM_fast32.argtypes     = TLBM_fast_argtypes
        TLBM_fast                          = TLBM_fast.TLBM_fast32

elif float_bits == 64:
     if RKTLBM:
        TLBM_fast.RKTLBM_strmac64.argtypes = RKTLBM_fast_strmac_argtypes
        RKTLBM_strmac                      = TLBM_fast.RKTLBM_strmac64
        TLBM_fast.RKTLBM_rest64.argtypes   = RKTLBM_fast_rest_argtypes
        RKTLBM_rest                        = TLBM_fast.RKTLBM_rest64
     else:
        TLBM_fast.TLBM_fast64.argtypes     = TLBM_fast_argtypes
        TLBM_fast                          = TLBM_fast.TLBM_fast64

#===============================================================================
# LBM parameter initialization for 2D (D2Q9) or 3D (D3Q19)
#===============================================================================
# Start of LBM parameter initialization
#-------------------------------------------------------------------------------


if ny == 1:  D = 2       # 2D case
else:        D = 3       # 3D case

D   = np_int(D)

R   = 1                  # Gas constant
C_s = 1/np.sqrt(3)       # Acoustic wave speed

DR_2 = np_float(D*R/2)   # Constant to convert between eta and temperature

# Velocity vectors for 2D and 3D, and weights

if D == 2:
    w0 = 4/9
    w1 = 1/9
    w2 = 1/36
    w  =  np.array((w0,w1,w1,w1,w1,w2,w2,w2,w2),dtype=np_float)
        
    B0 = -4/27
    B1 =  2/27
    B2 =  5/108
    B  =  np.array((B0,B1,B1,B1,B1,B2,B2,B2,B2),dtype=np_float)

    na = 9      # D2Q9
#   Symbolic names -> directions. ie. O0 = zero vector, R0=right, L0=left, U0=up, D0=down, RU=right-up, LD=left-down etc.
    O0       =  0               ; R0 =1 ; D0 =2 ; L0 =3 ; U0 =4 ; RD =5 ; LD =6 ; LU =7 ; RU =8
    c        = np.array((( 0, 0),( 1, 0),( 0, 1),(-1, 0),( 0,-1),( 1, 1),(-1, 1),(-1,-1),( 1,-1)),dtype=np_int)

#                           0       1       2       3       4       5       6       7       8
    ai       = np.array((   0    ,  3   ,   4   ,   1   ,   2   ,   7   ,   8   ,   5   ,   6  ),dtype=np_int)
elif D == 3:
    w0 = 1/3
    w1 = 1/18
    w2 = 1/36
    w  =  np.array((w0,w1,w1,w1,w1,w1,w1,w2,w2,w2,w2,w2,w2,w2,w2,w2,w2,w2,w2),dtype=np_float)
    na = 19     # D3Q19
#   Symbolic names -> directions. ie. O0 = zero vector, R0=right, L0=left, U0=up, D0=down, RU=right-up, LD=left-down etc.
    O0       =              0       ;  R0 =1    ;  L0 =2    ;  F0 =3    ;  B0 =4    ;  D0 =5    ;  U0 =6    ;  RF =7    ;  LF =8    ;  RB =9    ;  LB =10   ;   RD =11  ;  LD =12    ;  RU =13   ;  LU =14   ;  DF =15   ;  DB =16   ;  UF =17 ;  UB = 18
    c        = np.array([ [ 0, 0, 0], [ 1, 0, 0], [-1, 0, 0], [ 0, 1, 0], [ 0,-1, 0], [ 0, 0, 1], [ 0, 0,-1], [ 1, 1, 0], [-1, 1, 0], [ 1,-1, 0], [-1,-1, 0], [ 1, 0, 1], [ -1, 0, 1], [ 1, 0,-1], [-1, 0,-1], [ 0, 1, 1], [ 0,-1, 1], [ 0, 1,-1], [0, -1, -1] ],dtype=np_int)
#                           0           1           2           3           4           5           6           7           8           9          10          11           12          13          14          15          16          17          18
    ai       = np.array([   0       ,   2       ,   1       ,   4       ,   3       ,   6       ,   5       ,  10       ,   9       ,   8       ,   7       ,  14       ,   13       ,  12       ,  11       ,  18       ,  17       ,  16       ,  15      ],dtype=np_int)
    
    #   3D coefficients for B from: Liu, H., Valocchi, A. J., and Kang, Q. (2012)
    #   Three-Dimensional Lattice Boltzmann Model for Immiscible Two-Phase Flow Simulations
    #   Physical Review E. 71: 046309, https://doi.org/10.1103/PhysRevE.85.046309.

    B0 = -1/3
    B1 = 1/18
    B2 = 1/36
    B  =  np.array((B0,B1,B1,B1,B1,B1,B1,B2,B2,B2,B2,B2,B2,B2,B2,B2,B2,B2,B2),dtype=np_float)

# MRT matrices M and vector S in 2D and 3D

if D == 2:
    M = np.array([[ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0],
                  [-4.0,-1.0,-1.0,-1.0,-1.0, 2.0, 2.0, 2.0, 2.0],
                  [ 4.0,-2.0,-2.0,-2.0,-2.0, 1.0, 1.0, 1.0, 1.0],
                  [ 0.0, 1.0, 0.0,-1.0, 0.0, 1.0,-1.0,-1.0, 1.0],
                  [ 0.0,-2.0, 0.0, 2.0, 0.0, 1.0,-1.0,-1.0, 1.0],
                  [ 0.0, 0.0, 1.0, 0.0,-1.0, 1.0, 1.0,-1.0,-1.0],
                  [ 0.0, 0.0,-2.0, 0.0, 2.0, 1.0, 1.0,-1.0,-1.0],
                  [ 0.0, 1.0,-1.0, 1.0,-1.0, 0.0, 0.0, 0.0, 0.0],
                  [ 0.0, 0.0, 0.0, 0.0, 0.0, 1.0,-1.0, 1.0,-1.0]],dtype=np_float)
    Mi = inv(M)

    tau_f = 1        # Just set to something here. Later it will be set
                     # such that tau_f and hence viscosity are variable
    omega = 1/tau_f
    s0 = 1      ; s3 = s0 ; s5 = s0 # rho, j_x, j_y
    s1 = 1.4    ; s2 = s1           # e and eta
    s4 = 1.2    ; s6 = s4           # q_x and q_y
    s7 = omega  ; s8 = s7           # p_xx and p_xy
    S_  = np.array([s0,s1,s2,s3,s4,s5,s6,s7,s8],dtype=np_float)

#-------------------------------------------------------------------------------
# Constants for M & S of the 2D MRT algorithm are taken from the paper:
#
# E. Aslan, I. Taymaz, and A. C. Benim, 2014, Investigation of the Lattice
# Boltzmann SRT and MRT Stability for Lid Driven Cavity Flow, Int. J. of
# Materials, Mechanics and Manufacturing, Vol. 2, No. 4, November 2014
#
# Parameters s0 -> s6 are "free parameters" and can be changed to optimize
# stability, whereas s7 and s8 must equal to 1/tau_f. As tau_f may be
# spatially variable, this means we first initialize S as 1d vector, then
# later we set s7 and s8 to the calculated value of 1/tau_f(x,t) to get the
# correct viscosity.
#-------------------------------------------------------------------------------

elif D == 3:
    M = np.array([[  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0],  # 0
                  [-30.0,-11.0,-11.0,-11.0,-11.0,-11.0,-11.0,  8.0,  8.0,  8.0,  8.0,  8.0,  8.0,  8.0,  8.0,  8.0,  8.0,  8.0,  8.0],  # 1
                  [ 12.0, -4.0, -4.0, -4.0, -4.0, -4.0, -4.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0],  # 2
                  [  0.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0],  # 3
                  [  0.0, -4.0,  4.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0],  # 4
                  [  0.0,  0.0,  0.0,  1.0, -1.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0],  # 5
                  [  0.0,  0.0,  0.0, -4.0,  4.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0],  # 6
                  [  0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0, -1.0],  # 7
                  [  0.0,  0.0,  0.0,  0.0,  0.0, -4.0,  4.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0, -1.0],  # 8
                  [  0.0,  2.0,  2.0, -1.0, -1.0, -1.0, -1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0, -2.0, -2.0, -2.0, -2.0],  # 9
                  [  0.0, -4.0, -4.0,  2.0,  2.0,  2.0,  2.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0, -2.0, -2.0, -2.0, -2.0],  #10
                  [  0.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0,  1.0,  1.0,  1.0,  1.0, -1.0, -1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0],  #11
                  [  0.0,  0.0,  0.0, -2.0, -2.0,  2.0,  2.0,  1.0,  1.0,  1.0,  1.0, -1.0, -1.0, -1.0, -1.0,  0.0,  0.0,  0.0,  0.0],  #12
                  [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0, -1.0,  1.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0],  #13
                  [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0, -1.0,  1.0],  #14
                  [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0, -1.0,  1.0,  0.0,  0.0,  0.0,  0.0],  #15
                  [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0, -1.0,  1.0, -1.0,  1.0,  0.0,  0.0,  0.0,  0.0],  #16
                  [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0, -1.0, -1.0,  1.0,  1.0,  0.0,  0.0,  0.0,  0.0,  1.0, -1.0,  1.0, -1.0],  #17
                  [  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  0.0,  1.0,  1.0, -1.0, -1.0, -1.0, -1.0,  1.0,  1.0]],dtype=np_float) #18
    Mi = inv(M)

    tau_f = 1        # Just set to something here. Later it will be set
                     # such that tau_f and hence viscosity are variable
    s1 = 1.0; s4 = s1; s6 = s1;  s8 = s1
    s2 = 1.19
    s3 = 1.4; s11=s3; s13=s3
    s5 = 1.2; s7 = s5; s9 = s5
    s17 = 1.98; s18=s17; s19=s7
    s10 = 1/tau_f  ; s12 = s10; s14=s10; s15=s10; s16=s10
    S_  = np.array([s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15,s16,s17,s18,s19],dtype=np_float)

#-------------------------------------------------------------------------------
# Constants for M & S of the 3D MRT algorithm are taken from the paper:
#
# Li, Zheng, Mo Yang, and Yuwen Zhang. "Lattice Boltzmann method simulation of
# 3-D natural convection with double MRT model."
# International Journal of Heat and Mass Transfer 94 (2016): 222-238.
#
#-------------------------------------------------------------------------------

Mi_S = np.zeros((na,na),dtype=np_float) # oooOOO Should be done during the simulation as tau_f and hence omega changes
for a in range(na):
    for b in range(na):
        Mi_S[a,b] += Mi[a,b] * S_[b]

#-------------------------------------------------------------------------------
# End of LBM parameter initialization
#===============================================================================
# Run initialization from the input specifications
#===============================================================================
# Start of run specifications
#-------------------------------------------------------------------------------

if rank == 0: print("Parameter initialization ...",flush=True)

# Initialize latex if necessary

if latex:
    matplotlib.rcParams['text.usetex'] = True
else:
    matplotlib.rcParams['text.usetex'] = False

#-------------------------------------------------------------------------------
# Thermal BC, lattice size, kinematic viscosity, geometry etc
#-------------------------------------------------------------------------------

tiny = np.finfo(np_float).tiny ; zero = tiny  # The smallest float for dtype = np_float

# Temperature at boundaries and constants for temperature dependent viscosity, kinematic viscosity

T_upper  = -dT_p/2      # Upper temperature
T_lower  =  dT_p/2      # Lower temperature

# Initialize for the desired geometry

R_p_input = R_p
L_p_input = L_p
if geometry == "annular" or (geometry == "spherical" and D == 2):
    L_l,nz,R_p,L_p,radii,thetas,nu_factor,dy_ann,dx_ann = get_annular_geometry(R_p,L_p,theta1,theta2,nx,rank)
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

dx_p     = L_p/L_l      # Size of lattice,      m

# Initialize the boundary indices and lattice size
# Note: it is lattice spacing at the surface for annular geometry, nu_factor takes care of different spacing in r

z_upper  = 0            # Upper boundary where the temperature BC is set
z_lower  = nz-1         # Lower boundary where the temperature BC is set

# Select default plane and coordinate for 2D plot_field function

for i in range(len(coordinates)):
    coordinate = int(coordinates[i])
    if   plane == "xz" and (coordinate < 0 or coordinate >= ny):
        print("Warning: plane=%s, coordinate=%d not in vaid range of 0 to %d, being set to %d" %(plane,coordinate,ny-1,ny//2),flush=True)
        coordinate = ny//2
        coordinates[i] = str(coordinate)
    elif plane == "xy" and (coordinate < 0 or coordinate >= nz):
        print("Warning: plane=%s, coordinate=%d not in vaid range of 0 to %d, being set to %d" %(plane,coordinate,nz-1,nz//2),flush=True)
        coordinate = nz//2
        coordinates[i] = str(coordinate)
    elif plane == "yz" and (coordinate < 0 or coordinate >= nx):
        print("Warning: plane=%s, coordinate=%d not in vaid range of 0 to %d, being set to %d" %(plane,coordinate,nx-1,nx//2),flush=True)
        coordinate = nx//2
        coordinates[i] = str(coordinate)

#-------------------------------------------------------------------------------
# Calculate the lattice viscosity, thermal diffusivity, gravity and
# physical time step, and calculate Ra_p and Pr_p
#-------------------------------------------------------------------------------

# Convert to secs if timescale != unspecified

time_per_timestep_input = time_per_timestep
time_per_timestep_int, timescale_time_per_timestep, error = timevalue("time_per_timestep",time_per_timestep,error,rank)
time_per_timestep = timevalue_to_secs(time_per_timestep_int,timescale_time_per_timestep,rank)

eta_p_input = eta_p; nu_l_input = nu_l ; g_l_input = g_l
Ra_p, Pr_p, eta_p, nu_p, ka_p, Ra_l, Pr_l, g_l, nu_l, ka_l, dx_p, dt_p, dt_dt0 = get_lattice_viscosity(state, power_correction,Ra,Pr,g_p,beta_p,dT_p,L_p,eta_p,rho_p,ka_p,Tau_y0_p,D,g_l,L_l,nu_l,nu_l_min,time_per_timestep,gamma_T,rheology,rank)

eta_p_input = np_float(eta_p_input)
# Initialize the viscosities for distributions f and g

nu_f0        = nu_l
nu_g0        = ka_l

if RKTLBM:
   nu_g0    = nu_f0/Pr
   # ka_l     = nu_g0
   nu_f0_r  = nu_f0   * nu_r_factor # Red  viscosity
   nu_f0_b  = nu_f0                 # Blue viscosity
   tau_f0_r = nu_f0_r * 3 + 0.5     # Default relaxation time for red
   tau_f0_b = nu_f0_b * 3 + 0.5     # Default relaxation time for blue

if RKTLBM:
   if D == 2: 
      beta_g = np_float(-3 * beta_p * g_l)                     # Boussinesq factor
   else:
      beta_g = np_float(-3 * beta_p * g_l)                     # Boussinesq factor (1/cs², same as D=2 and non-RKTLBM)
   if rho1 == rho2 and rho1_p != rho2_p:
      # Lattice densities left at defaults would silently disable compositional
      # buoyancy (g_comp = 0). Auto-derive them from the physical densities.
      rho1 = np_float(rho1_p/rho2_p)
      rho2 = np_float(1.0)
      if rank == 0:
          print("NOTE: rho1 == rho2 with rho1_p != rho2_p would give g_comp = 0.",flush=True)
          print("      Auto-setting lattice rho1 = rho1_p/rho2_p = %.6g, rho2 = 1.0" %(rho1),flush=True)
          print("      (If RHO_r/RHO_b are read from .npy files, they must use these values.)",flush=True)
   if rho1 != rho2:
      g_comp = g_l * (rho1_p-rho2_p)*rho2/(rho2_p*(rho1-rho2)) # Compositional gravity
   else: 
      g_comp = 0 
   g_comp   = np_float(g_comp)
   beta_g  *= beta_g_factor
   g_comp  *= g_comp_factor

else: beta_g   = np_float(-3 * beta_p * g_l)

# Initialize the relaxation times for distributions f and g

tau_f0    = nu_f0   * 3 + 0.5              # Initialize relaxation time for f
tau_g0    = nu_g0   * 3 + 0.5              # Initialize relaxation time for g
tau_h0    = nu_g0   * Pr/Peclet * 3 + 0.5  # Initialize relaxation time for h

# Ensure all values are at the correct precision

nu_f0  = np_float(nu_f0)    ; tau_f0   = np_float(tau_f0)
nu_g0  = np_float(nu_g0)    ; tau_g0   = np_float(tau_g0)
tau_h0 = np_float(tau_h0)

if RKTLBM:
   nu_f0_r   = np_float(nu_f0_r)  ; tau_f0_r = np_float(tau_f0_r)
   nu_f0_b   = np_float(nu_f0_b)  ; tau_f0_b = np_float(tau_f0_b)

nu_l_min  = np_float(nu_l_min)

# Physical viscosity clamps (Pa.s) -> lattice kinematic viscosity.
# Inverse of the forward conversion used in the nu_kinematic/dynamic plots:
#     nu_p = (Pr_p/Pr_l) * nu_l * dx_p^2 / dt_p    (state=solid)
#     nu_p = sqrt(Ra_l/Ra_p) * nu_l * dx_p^2 / dt_p (else)
# nu_p = eta_p/rho_p
if state == "solid":
    _nu_unit = (Pr_p/Pr_l) * dx_p**2 / dt_p
else:
    _nu_unit = np.sqrt(Ra_l/Ra_p) * dx_p**2 / dt_p
nu_l_min_cap = np_float((eta_min_p / rho_p) / _nu_unit)
nu_l_max_cap = np_float((eta_max_p / rho_p) / _nu_unit)

# theta < 0 selects neutral (zero-gradient) wall ghosts for rho_r/rho_b:
# the ghost keeps the adjacent cell's red/blue split instead of imposing a
# fixed wetting fraction, which otherwise pumps color at free surfaces.
BC_r_wetting_u = np_float(-1.0 if theta_u < 0 else (1+np.cos(np.pi*theta_u/180))/2)
BC_b_wetting_u = np_float(-1.0 if theta_u < 0 else (1-np.cos(np.pi*theta_u/180))/2)
BC_r_wetting_d = np_float(-1.0 if theta_d < 0 else (1+np.cos(np.pi*theta_d/180))/2)
BC_b_wetting_d = np_float(-1.0 if theta_d < 0 else (1-np.cos(np.pi*theta_d/180))/2)

#==========================================================================
# Print the run specifications
#==========================================================================

nt_timesteps = timevalue_to_timesteps(nt,timescale_nt,dt_p,rank)

if rank == 0:
    print("==================================================================",flush=True)
    print("Run title and description on plots",flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("run_title                     = ",run_title,flush=True)
    print("description                   = ",description,flush=True)
    print("description_color             = ",description_color,flush=True)
    print("description_box               = ",description_box,flush=True)
    print("description_location          = ",description_location,flush=True)
    run_title   = str2str(run_title)               # Change underscores to spaces except those preceded by an underscore
    description = str2str(description)             # Change underscores to spaces except those preceded by an underscore
    print("------------------------------------------------------------------",flush=True)
    print("Model size, timesteps and geometry",flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("nx                           = ",nx,flush=True)
    print("ny                           = ",ny,flush=True)
    print("nz                           = ",nz,flush=True)
    print("L_l                          = %8.5g" %(L_l),flush=True)
    print("nt                           = %d %s = %d timesteps" %(nt,timescale_nt,nt_timesteps),flush=True)
    print("geometry                     = ",geometry,flush=True)
    if type(R_p_input) == str:
        print("R_p                          = ",R_p_input,flush=True)
    print("R_p         (m)              = %8.3g" %(R_p),flush=True)
    if type(L_p_input) == str:
        print("L_p                          = ",L_p_input,flush=True)
    print("L_p         (m)              = %8.3g" %(L_p),flush=True)
    print("theta1      (degrees)        = %8.3g" %(theta1),flush=True)
    print("theta2      (degrees)        = %8.3g" %(theta2),flush=True)
    print("n_TBL                        = %d" %(n_TBL),flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("Physical specifications",flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("state                        = ",state,flush=True)
    print("power_correction             = ",power_correction,flush=True)
    if Ra == "unspecified":
        print("Ra                           = ",Ra,flush=True)
    else:
        print("Ra                           = %8.3g" %(Ra),flush=True)
    print("Pr                           = %8.3g" %(Pr),flush=True)
    print("L_p         (m)              = %8.3g" %(L_p),flush=True)
    print("rho_p       (kg/m^3)         = %8.3g" %(rho_p),flush=True)
    print("Density profile file         = ", rho_z_file,flush=True)
    print("Temperature profile file     = ", T_z_file,flush=True)
    print("Pressure profile file        = ", P_z_file,flush=True)
    print("Lithology profile file       = ", L_z_file,flush=True)
    if eta_p_input == "unspecified":
        print("eta_p                        = ",eta_p_input,flush=True)
    else:
        print("eta_p       (kg/(m.s))       = %8.3g" %(eta_p_input),flush=True)
    print("eta_min_p   (Pa.s)           = %8.3g" %(eta_min_p),flush=True)
    print("eta_max_p   (Pa.s)           = %8.3g" %(eta_max_p),flush=True)
    print("act_E_1     (J/mol)          = %8.3g" %(act_E_1),flush=True)
    print("act_E_2     (J/mol)          = %8.3g" %(act_E_2),flush=True)
    print("act_E_3     (J/mol)          = %8.3g" %(act_E_3),flush=True)
    print("act_V_1     (m^3/mol)        = %8.3g" %(act_V_1),flush=True)
    print("act_V_2     (m^3/mol)        = %8.3g" %(act_V_2),flush=True)
    print("act_V_3     (m^3/mol)        = %8.3g" %(act_V_3),flush=True)
    print("Delta_eta_1                  = %8.3g" %(Delta_eta_1),flush=True)
    print("Delta_eta_2                  = %8.3g" %(Delta_eta_2),flush=True)
    print("Delta_eta_3                  = %8.3g" %(Delta_eta_3),flush=True)
    print("P_scale_1   (Pa)             = %8.3g" %(P_scale_1),flush=True)
    print("P_scale_2   (Pa)             = %8.3g" %(P_scale_2),flush=True)
    print("P_scale_3   (Pa)             = %8.3g" %(P_scale_3),flush=True)
    print("ka_p        (m^2/s)          = %8.3g" %(ka_p),flush=True)
    print("g_p         (m/s^2)          = ",g_p,flush=True)
    print("beta_p      (1/K)            = %8.3g" %(beta_p),flush=True)
    print("T0_p        (C)              = ",T0_p,flush=True)
    print("dT_p        (C)              = ",dT_p,flush=True)
    print("dT_p_random (C)              = ",dT_p_random,flush=True)
    print("T_p_random  (C)              = ",T_p_random,flush=True)
    print("initialize_ncells            = ",initialize_ncells,flush=True)
    print("cell_amplitude               = ",cell_amplitude,flush=True)
    print("amp_sinusoid                 = ",amp_sinusoid,flush=True)
    print("x_periods                    = ",x_periods,flush=True)
    print("y_periods                    = ",y_periods,flush=True)
    print("H                            = ",H,flush=True)
    print("Cp                           = ",Cp,flush=True)
    print("temperature                  = ",temperature,flush=True)
    print("gamma_T     (1/K)            = ",gamma_T,flush=True)
    print("rheology                     = ",rheology,flush=True)
    print("nu_formula                   = ",nu_formula,flush=True)
    print("Tau_y0_p    (Pa=kg/(m.s^2))  = %8.3g" %(Tau_y0_p),flush=True)
    print("Tau_yi_p    (Pa/m)           = %8.3g" %(Tau_yi_p),flush=True)
    print("nu_star                      = %8.3g" %(nu_star),flush=True)
    print("N           (exponent)       = ",N,flush=True)
    print("RKTLBM                       = ",RKTLBM,flush=True)
    print("collision                    = ",collision,flush=True)
    print("RKTLBM red density file      = ",rho_r_file,flush=True)
    print("RKTLBM blue density file     = ",rho_b_file,flush=True)
    print("RKTLBM R relative density (rho1) = ",rho1,flush=True)
    print("RKTLBM B relative density (rho2) = ",rho2,flush=True)
    print("RKTLBM red physical density  = ",rho1_p,flush=True)
    print("RKTLBM blue physical density = ",rho2_p,flush=True)
    print("RKTLBM R/B density (rho12)   = ",rho12,flush=True)
    print("A                            = ",A,flush=True)
    print("interface_thickness          = ",interface_thickness,flush=True)
    print("nu_r_factor                  = ",nu_r_factor,flush=True)
    print("theta_u                      = ",theta_u)
    print("theta_d                      = ",theta_d)
    print("basalt_depth_p  (m)          = %8.3g" %(basalt_depth_p),flush=True)
    print("BC_r_wetting_u               = ",BC_r_wetting_u,flush=True)
    print("BC_r_wetting_d               = ",BC_r_wetting_d,flush=True)
    print("BC_b_wetting_u               = ",BC_b_wetting_u,flush=True)
    print("BC_b_wetting_d               = ",BC_b_wetting_d,flush=True)
    if RKTLBM:
       print("g_comp_factor                = %8.3g" %(g_comp_factor),flush=True)
       print("beta_g_factor                = %8.3g" %(beta_g_factor),flush=True)
       print("g_comp                       = %8.3g" %(g_comp),flush=True)
       print("beta_g                       = %8.3g" %(beta_g),flush=True)
       print("melt_time                    = ",melt_time,flush=True)
       print("melt_position                = ",melt_position,flush=True)
       print("melt_extraction              = ",melt_extraction,flush=True)
       print("melt_test_mode               = ",melt_test_mode,flush=True)
    if rheology == "linear":
        print("   Tau_y0_p, Tau_yi_p and N are not applicable for rheology = linear",flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("Boundary condition specifications",flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("BC_sides                     = ",BC_sides,flush=True)
    print("BC_upper                     = ",BC_upper,flush=True)
    print("BC_lower                     = ",BC_lower,flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("Physical timestep specifications",flush=True)
    print("------------------------------------------------------------------",flush=True)
    if time_per_timestep_input == "unspecified":
        print("time_per_timestep            = ",time_per_timestep_input,flush=True)
    else:
        print("time_per_timestep            = ",time_per_timestep_input,flush=True)
        print("                             = %d %s" %(time_per_timestep_int,timescale_time_per_timestep),flush=True)
    if nu_l_input == "unspecified":
        print("nu_l                         = ",nu_l_input,flush=True)
    else:
        print("nu_l                         = %8.3g" %(nu_l_input),flush=True)
    if g_l_input == "unspecified":
        print("g_l                          = ",g_l_input,flush=True)
    else:
        print("g_l                          = %8.3g" %(g_l_input),flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("Mach specifications to maximize speed & control accuracy/stability",flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("Inc_vary_dt                  = %d %s" %(Inc_vary_dt  ,timescale_Inc_vary_dt  ),flush=True)
    print("start_vary_dt                = %d %s" %(start_vary_dt,timescale_start_vary_dt),flush=True)
    print("    factor_increase_dt       = ",factor_increase_dt,flush=True)
    print("    factor_decrease_dt       = ",factor_decrease_dt,flush=True)
    print("    factor_max_vary_dt       = ",factor_max_vary_dt,flush=True)
    print("    max_mach_limit           = ",max_mach_limit,flush=True)
    print("    min_mach_limit           = ",min_mach_limit,flush=True)
    print("nu_l_min                     = ",nu_l_min,flush=True)
    print("mach_limit_to_exit           = ",mach_limit_to_exit,flush=True)
    print("==================================================================",flush=True)
    print("Physical+lattice parameters, lattice spacing, & physical timestep",flush=True)
    print("==================================================================",flush=True)
    print("eta_p      (m^2/s)           = %8.3g" %(eta_p),flush=True)
    print("nu_p       (m^2/s)           = %8.3g" %(nu_p),flush=True)
    print("ka_p       (m^2/s)           = %8.3g" %(ka_p),flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("nu_l                         = %8.3g" %(nu_l),flush=True)
    print("ka_l                         = %8.3g" %(ka_l),flush=True)
    print("g_l                          = %8.3g" %(g_l),flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("dx_p = L_p/L_l (m)           = %8.3g" %(dx_p),flush=True)
    print("dt_p           (s)           = %8.3g" %(dt_p),flush=True)
    print("                             = %s" %(time_to_physical_time_text(dt_p)),flush=True)
    print("==================================================================",flush=True)
    print("Checks of TLBM values being used",flush=True)
    print("==================================================================",flush=True)
    print("nu_f = %8.3g , tau_f = %8.3g" %(nu_f0,tau_f0),flush=True)
    print("nu_g = %8.3g , tau_g = %8.3g" %(nu_g0,tau_g0),flush=True)
    print("nu_h = %8.3g , tau_h = %8.3g" %(nu_g0,tau_h0),flush=True)
    if RKTLBM:    
       print("nu_f_r = %8.3g" %(nu_f0_r),flush=True)
       print("nu_f_b = %8.3g" %(nu_f0_b),flush=True)
    Ra_check  = g_l*beta_p*dT_p*L_l**3/(nu_f0*nu_g0)
    print("Ra_l                         = %8.3g" %(Ra_check),flush=True)
    print("==================================================================",flush=True)
    print("Calculation and output accuracy specifications",flush=True)
    print("==================================================================",flush=True)
    print("shifted_distribution         = ",shifted_distribution,flush=True)
    print("float_bits                   = ",float_bits,flush=True)
    print("edge_bits                    = ",edge_bits,flush=True)
    print("data_bits                    = ",data_bits,flush=True)
    print("==================================================================",flush=True)
    print("Output specifications",flush=True)
    print("==================================================================",flush=True)
    print("output_plots                 = ",output_plots,flush=True)
    print("output_data                  = ",output_data,flush=True)
    print("output_timing                = ",output_timing,flush=True)
    print("output_units                 = ",output_units,flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("Output increment specifications",flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("Inc_plots                    = %d %s" %(Inc_plots ,timescale_Inc_plots) ,flush=True)
    print("Inc_data                     = %d %s" %(Inc_data  ,timescale_Inc_data)  ,flush=True)
    print("Inc_mach                     = %d %s" %(Inc_mach  ,timescale_Inc_mach)  ,flush=True)
    print("Inc_timing                   = %d %s" %(Inc_timing,timescale_Inc_timing),flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("Output plane specifications",flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("plane                        = ",plane,flush=True)
    print("coordinates                  = ",coordinates,flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("Output options",flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("plot_size                    = ",plot_size,flush=True)
    print("plot_size_factor             = ",plot_size_factor,flush=True)
    print("plot_colorbar                = ",plot_colorbar,flush=True)
    print("colormap                     = ",colormap,flush=True)
    print("latex                        = ",latex,flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("Output plot specifications",flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("plot_rho                     = ",plot_rho,flush=True)
    print("plot_T                       = ",plot_T,flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("plot_u                       = ",plot_u,flush=True)
    print("plot_nu_total                = ",plot_nu_total,flush=True)
    print("plot_nu_kinematic            = ",plot_nu_kinematic,flush=True)
    print("plot_nu_dynamic              = ",plot_nu_dynamic,flush=True)
    print("plot_Ra                      = ",plot_Ra,flush=True)
    print("plot_Tau                     = ",plot_Tau,flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("plot_vectors                 = ",plot_vectors,flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("Archean2 melt extract plots specifications",flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("plot_P                       = ",plot_P,flush=True)
    print("plot_L                       = ",plot_L,flush=True)
    print("plot_rho_eff                 = ",plot_rho_eff,flush=True)
    print("plot_rho_s                   = ",plot_rho_s,flush=True)
    print("plot_M_pd                    = ",plot_M_pd,flush=True)
    print("plot_M_bs                    = ",plot_M_bs,flush=True)
    print("plot_M_LP                    = ",plot_M_LP,flush=True)
    print("plot_M_MP                    = ",plot_M_MP,flush=True)
    print("plot_M_HP                    = ",plot_M_HP,flush=True)
    print("plot_ce_pd                   = ",plot_ce_pd,flush=True)
    print("plot_ce_bs                   = ",plot_ce_bs,flush=True)
    print("plot_ce_LP                   = ",plot_ce_LP,flush=True)
    print("plot_ce_MP                   = ",plot_ce_MP,flush=True)
    print("plot_ce_HP                   = ",plot_ce_HP,flush=True)
    print("plot_T_sol_pd                = ",plot_T_sol_pd,flush=True)
    print("plot_T_sol_bs                = ",plot_T_sol_bs,flush=True)
    print("plot_T_liq_pd                = ",plot_T_liq_pd,flush=True)
    print("plot_T_liq_bs                = ",plot_T_liq_bs,flush=True)
    print("plot_DI                      = ",plot_DI,flush=True)
    print("plot_rho_r                   = ",plot_rho_r,flush=True)
    print("plot_rho_b                   = ",plot_rho_b,flush=True)
    print("plot_FF_mag                  = ",plot_FF_mag,flush=True)
    print("plot_FF_x                    = ",plot_FF_x,flush=True)
    print("plot_FF_z                    = ",plot_FF_z,flush=True)
    print("plot_FF_y                    = ",plot_FF_y,flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("Output scale specifications",flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("vector_overlay               = ",vector_overlay,flush=True)
    print("vector_color                 = ",vector_color,flush=True)
    print("vector_interval              = ",vector_interval,flush=True)
    print("vector_length                = ",vector_length,flush=True)
    print("vector_width                 = ",vector_width,flush=True)
    print("annulus_interval             = ",annulus_interval,flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("clip_rho                     = ",clip_rho,flush=True)
    print("gpow_rho                     = ",gpow_rho,flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("clip_T                       = ",clip_T,flush=True)
    print("gpow_T                       = ",gpow_T,flush=True)
    print("fix_scale_T                  = ",fix_scale_T,flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("clip_u                       = ",clip_u,flush=True)
    print("gpow_u                       = ",gpow_u,flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("clip_nu_total                = ",clip_nu_total,flush=True)
    print("gpow_nu_total                = ",gpow_nu_total,flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("clip_Tau                     = ",clip_Tau,flush=True)
    print("gpow_Tau                     = ",gpow_Tau,flush=True)
    print("==================================================================",flush=True)
    print("Performance related specifications",flush=True)
    print("==================================================================",flush=True)
    print("code                         = ",code,flush=True)
    print("overlap                      = ",overlap,flush=True)
    print("domain decomposition         = ",domain_decomposition,flush=True)
    print("nxi_min                      = ",nxi_min,flush=True)
    print("nxi_max                      = ",nxi_max,flush=True)
    print("nyi_min                      = ",nyi_min,flush=True)
    print("nyi_max                      = ",nyi_max,flush=True)
    print("nzi_min                      = ",nzi_min,flush=True)
    print("nzi_max                      = ",nzi_max,flush=True)

#   Print out explanation of unit conversions

    print("==================================================================",flush=True)
    if state == "solid":
        print("Unit conversions for state = solid => Pr_l << Pr_p & Ra_l = Ra_p",flush=True)
    else:
        print("Unit conversions for state = magma => Ra_l << Ra_p & Pr_l = Pr_p",flush=True)
    print("==================================================================",flush=True)
    if state == "solid":
        c0 = Pr_p/Pr_l
        c1 = 1
    else:
        c0 = 1
        c1 = np.sqrt(Ra_l/Ra_p)
    print("c0 = Pr_p/Pr_l       = %8.3g" %(c0),flush=True)
    print("c1 = sqrt(Ra_l/Ra_p) = %8.3g" %(c1),flush=True)
    print("dx_p                 = %8.3g m" %(dx_p),flush=True)
    print("dt_p                 = %8.3g s" %(dt_p),flush=True)
    if rheology == "archean" or rheology == "archean2":
        conversion_Tau_p_2_l  = (1/rho_p)
        conversion_Tau_p_2_l *= (g_l/g_p) / dx_p**(D-1)
        print("------------------------------------------------------------------",flush=True)
        print("Stress units = F/A = m*a/dx**(D-1) where D = number of dimensions",flush=True)
        print("    stress    = (rho*dx**D)*dx/dt**2/dx**(D-1) = rho*dx**2/dt**2",flush=True)
        print("    rho_p                = %8.3g" %(rho_p),flush=True)
        print("    Tau_p_2_l = %8.3g = %8.3g = dt_p**2/(c0*rho_p*dx_p**2)" %(conversion_Tau_p_2_l,dt_p**2/(c0*rho_p*dx_p**2)),flush=True)
        print("    Tau_y0_p  = %8.3g" %(Tau_y0_p),flush=True)
        print("    Tau_y0_l  = %8.3g = Tau_p_2_l * Tau_y0_p" %(conversion_Tau_p_2_l * Tau_y0_p),flush=True)
        print("------------------------------------------------------------------",flush=True)
    print("g_p  = %8.3g = c0*g_l*dx_p/dt_p**2     where g_l  = %8.3g" %(c0*g_l*dx_p/dt_p**2,g_l),flush=True)
    print("nu_p = %8.3g = c0*c1*nu_l*dx_p**2/dt_p where nu_l = %8.3g" %(nu_l*c0*c1*dx_p**2/dt_p,nu_l),flush=True)
    print("ka_p = %8.3g = c1*ka_l*dx_p**2/dt_p    where ka_l = %8.3g" %(c1*ka_l*dx_p**2/dt_p,ka_l),flush=True)
    v_l       = 0.005
    v_p       = v_l*dx_p/dt_p
    v_p_cm_py = v_p * u_cm_py
    if state == "solid":
        print("v_p = v_l*dx_p/dt_p*secs/y*100 = %8.3g cm/y for v_l = %5.3f" %(v_p_cm_py,v_l),flush=True)
    else:
        print("v_p = v_l*dx_p/dt_p = %8.3g m/s for v_l = %5.3f" %(v_p,v_l),flush=True)

#--------------------------------------------------------------
# Convert from input time units to timesteps
#--------------------------------------------------------------

nt            = nt_timesteps
start_vary_dt = timevalue_to_timesteps(start_vary_dt,timescale_start_vary_dt,dt_p,rank)
Inc_vary_dt   = timevalue_to_timesteps(Inc_vary_dt  ,timescale_Inc_vary_dt  ,dt_p,rank)
Inc_plots     = timevalue_to_timesteps(Inc_plots    ,timescale_Inc_plots    ,dt_p,rank)
Inc_data      = timevalue_to_timesteps(Inc_data     ,timescale_Inc_data     ,dt_p,rank)
Inc_mach      = timevalue_to_timesteps(Inc_mach     ,timescale_Inc_mach     ,dt_p,rank)
Inc_timing    = timevalue_to_timesteps(Inc_timing   ,timescale_Inc_timing   ,dt_p,rank)

#==============================================================
# Initialization of domain decomposition
#==============================================================

# Get the domain decomposition

t0 = time.time()
if rank == 0:
    print("==================================================================",flush=True)
    print("Starting domain decomposition ...",end="",flush=True)

nr,nrx,nry,nrz,nxi,nyi,nzi,nx_i,ny_i,nz_i,mx,my,mz = my_domain_decomposition(size,rank,nx,ny,nz,nxi_min,nxi_max,nyi_min,nyi_max,nzi_min,nzi_max,domain_decomposition)

Dt = time.time() - t0
if rank == 0: print(" %4.2f secs" %(Dt),flush=True)

if rank == 0:
    mx0     = mx[0]       ; my0     = my[0]       ; mz0     = mz[0]
    rx_last = max(0,nrx-1) ; ry_last = max(0,nry-1) ; rz_last = max(0,nrz-1)
    mx_last = mx[rx_last] ; my_last = my[ry_last] ;mz_last = mz[rz_last] ;
    print("------------------------------------------------------------------",flush=True)
    print("Domain decomposition",flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("                    x       y       z",flush=True)
    print("Model size:   = %5d   %5d   %5d" %(nx,ny,nz),flush=True)
    print("------------------------------------------------------------------",flush=True)
    print("# domains     = %5d   %5d   %5d" %(nrx,nry,nrz),flush=True)
    print("Domain sizes  = %5d   %5d   %5d" %(mx0,my0,mz0),flush=True)
    if mx_last != mx0 or my_last != my0 or mz_last != mz0:
        print("Edge    sizes = %5d   %5d   %5d" %(mx_last,my_last,mz_last),flush=True)
    print("------------------------------------------------------------------",flush=True)
    cores_used = nrx*nry*nrz
    cores_idle = size - cores_used
    if cores_idle == 0:
        print("Cores = %d, all cores are being used" %(size),flush=True)
    else:
        print("Cores       = %d" %(size),flush=True)
        print("Used  cores = %d" %(cores_used),flush=True)
        print("Idle  cores = %d" %(cores_idle),flush=True)
    print("------------------------------------------------------------------",flush=True)

# End of initialization of MPI domains
#==============================================================
# Start of array allocation

if rank == 0: print("Allocating arrays ... ",end="",flush=True)

t0 = time.time()

rho       = np.ones((ny_i,nz_i,nx_i),dtype=np_float)     # Macroscopic density
eta       = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)    # Macroscopic energy density
u         = np.zeros((ny_i,nz_i,D,nx_i),dtype=np_float)  # Macroscopic velocity

if code == "c":
   if RKTLBM: 
      rho_r          = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)     # Red density
      rho_b          = np.ones((ny_i,nz_i,nx_i),dtype=np_float)      # Blue density
      rho_reff       = np.ones((ny_i,nz_i,nx_i),dtype=np_float)
      f_r            = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float)  # Mass density distribution
      f_r_last       = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float)  # Mass density distribution
      f_b            = np.ones((ny_i,nz_i,na,nx_i),dtype=np_float)   # Mass density distribution
      f_b_last       = np.ones((ny_i,nz_i,na,nx_i),dtype=np_float)
      g              = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float)  # Energy density distribution
      g_last         = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float)  # Energy density distribution
      FF             = np.zeros((D,ny_i,nz_i,nx_i),dtype=np_float)
      tau_g          = tau_g0 * np.ones((ny_i,nz_i,nx_i),dtype=np_float)
      tau_h          = tau_h0 * np.ones(nz_i,dtype=np_float)
      nu_f_b         = nu_f0_b * np.ones((nz_i),dtype=np_float)
      nu_f_r         = nu_f0_r * np.ones((nz_i),dtype=np_float)
      # Archean2 arrays
      h              = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float)  # Lithology density distribution
      h_last         = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float)  # Lithology density distribution
      hd             = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float)  # Advected depletion distribution (peridotite)
      hd_last        = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float)  # Advected depletion distribution (peridotite)
      hb             = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float)  # Advected depletion distribution (basalt)
      hb_last        = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float)  # Advected depletion distribution (basalt)
      L_ed                    = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      D_ed                    = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      D_eb                    = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      T_nd                    = np.ones((ny_i,nz_i,nx_i),dtype=np_float)
      z_nd                    = np.ones((ny_i,nz_i,nx_i),dtype=np_float) 
      M_ext_pd                = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      M_ext_bs                = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      M_ext_bs_LP_TTG         = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      M_ext_bs_MP_TTG         = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      M_ext_bs_HP_TTG         = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      M_pd                    = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      M_bs                    = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      cumulative_extracted_pd = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      cumulative_extracted_bs = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      cumulative_extracted_LP = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      cumulative_extracted_MP = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      cumulative_extracted_HP = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      rho_solid               = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      rho_eff                 = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      rho_M                   = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      P0                      = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      mask_pd                 = np.ones((ny_i,nz_i,nx_i),dtype=bool)
      mask_bs                 = np.ones((ny_i,nz_i,nx_i),dtype=bool)
      nu_total                = np.ones((ny_i,nz_i,nx_i),dtype=np_float)    # Total viscosity at last time step     
      Cex_BS                  = np.zeros((ny,nz,nx),dtype=np.float64)
      Cex_PD                  = np.zeros((ny,nz,nx),dtype=np.float64)
      Cex_LP                  = np.zeros((ny,nz,nx),dtype=np.float64)
      Cex_MP                  = np.zeros((ny,nz,nx),dtype=np.float64)
      Cex_HP                  = np.zeros((ny,nz,nx),dtype=np.float64)
   else:
      f         = np.ones((ny_i,nz_i,na,nx_i),dtype=np_float)  # Mass density distribution
      f_last    = np.ones((ny_i,nz_i,na,nx_i),dtype=np_float)  # Mass density distribution
      g         = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float) # Energy density distribution
      g_last    = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float) # Energy density distribution
      h         = np.ones((ny_i,nz_i,na,nx_i),dtype=np_float)  # Lithology density distribution
      h_last    = np.ones((ny_i,nz_i,na,nx_i),dtype=np_float)  # Lithology density distribution
      drho      = np.zeros(nx_i,dtype=np_float)                # Density perturbations
      f_eq      = np.zeros((na,nx_i),dtype=np_float)           # Equilibrium distribution f
      f_neq     = np.zeros((na,nx_i),dtype=np_float)           # Non-equilibrium distribution f
      g_neq     = np.zeros((na,nx_i),dtype=np_float)           # Non-equilibrium distribution g
      nu_f      = nu_f0    * np.ones(nz_i,dtype=np_float)         # Kinematic viscosity vs depth
      nu_g      = nu_g0    * np.ones(nz_i,dtype=np_float)         # Thermal diffusivity vs depth
      tau_g     = tau_g0   * np.ones(nz_i,dtype=np_float)
      tau_h     = tau_h0   * np.ones(nz_i,dtype=np_float)         # Relaxation time h vs depth
      # Archean2 arrays
      h_eq                    = np.zeros((na,nx_i),dtype=np_float)           # Equilibrium distribution h
      h_neq                   = np.zeros((na,nx_i),dtype=np_float)           # Non-equilibrium distribution h
      M_ext_pd                = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      M_ext_bs                = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      M_ext_bs_LP_TTG         = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      M_ext_bs_MP_TTG         = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      M_ext_bs_HP_TTG         = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      M_pd                    = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      M_bs                    = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      cumulative_extracted_pd = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      cumulative_extracted_bs = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      cumulative_extracted_LP = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      cumulative_extracted_MP = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      cumulative_extracted_HP = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      rho_solid               = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      rho_eff                 = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      rho_M                   = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      P0                      = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      L_ed                    = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
      T_nd                    = np.ones((ny_i,nz_i,nx_i),dtype=np_float)
      z_nd                    = np.ones((ny_i,nz_i,nx_i),dtype=np_float)
      mask_pd                 = np.ones((ny_i,nz_i,nx_i),dtype=bool)
      mask_bs                 = np.ones((ny_i,nz_i,nx_i),dtype=bool)
      nu_total   = np.ones((ny_i,nz_i,nx_i),dtype=np_float)    # Total viscosity at last time step
      
if code == "python":
    f        = np.ones((ny_i,nz_i,na,nx_i),dtype=np_float)  # Mass density distribution
    f_last   = np.ones((ny_i,nz_i,na,nx_i),dtype=np_float)  # Mass density distribution
    g        = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float) # Energy density distribution
    g_last   = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float) # Energy density distribution
    h        = np.ones((ny_i,nz_i,na,nx_i),dtype=np_float)  # Lithology density distribution
    h_last   = np.ones((ny_i,nz_i,na,nx_i),dtype=np_float)  # Lithology density distribution
    f_eq     = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float) # Equilibrium distribution f
    g_eq     = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float) # Non-equilibrium distribution g
    h_eq     = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float) # Equilibrium distribution h
    F_b      = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float) # Boussinesq term
    nu_f     = nu_f0   * np.ones((ny_i,nz_i,nx_i),dtype=np_float)  # Kinematic viscosity vs depth
    nu_g     = nu_g0   * np.ones((ny_i,nz_i,nx_i),dtype=np_float) # Thermal diffusivity vs depth
    if not (code == "c" and RKTLBM):
        nu_f_b = nu_f0_b * np.ones((ny_i,nz_i,nx_i),dtype=np_float)
    tau_f    = tau_f0  * np.ones((ny_i,nz_i,nx_i),dtype=np_float) # Relaxation time for f
    tau_g    = tau_g0  * np.ones((ny_i,nz_i,nx_i),dtype=np_float) # Relaxation time g vs depth
    nu_total = np.ones((ny_i,nz_i,nx_i),dtype=np_float)    # Total viscosity at last time step
    if RKTLBM:
       rho_r      = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)    # Red density
       rho_b      = np.ones((ny_i,nz_i,nx_i),dtype=np_float)     # Blue density
       rho_reff   = np.ones((ny_i,nz_i,nx_i),dtype=np_float) 
       f_r        = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float) # Mass density distribution
       f_r_last   = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float) # Mass density distribution
       f_b        = np.ones((ny_i,nz_i,na,nx_i),dtype=np_float)  # Mass density distribution
       f_b_last   = np.ones((ny_i,nz_i,na,nx_i),dtype=np_float)
       feq_r      = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float)
       feq_b      = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float)
       f_r_neq    = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float)
       f_b_neq    = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float)     
       F          = np.zeros((D,nx_i),dtype=np_float) 
       FF         = np.zeros((D,ny_i,nz_i,nx_i),dtype=np_float)
       fm         = np.zeros((nx_i),dtype=np_float)
       ttmp       = np.zeros((nx_i),dtype=np_float)
       cm         = np.zeros((na),dtype=np_float)
       cos_phi    = np.zeros((na,nx_i),dtype=np_float)
       cos_lambda = np.zeros((na,nx_i),dtype=np_float)       
       ff_sum     = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float)
       feq        = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float)
       tempp      = np.zeros((ny_i,nz_i,na,nx_i),dtype=np_float)


nu_f_last  = np.ones((ny_i,nz_i,nx_i),dtype=np_float)    # Kinematic viscosity at last time step
tau_f_last = np.ones((ny_i,nz_i,nx_i),dtype=np_float)    # Relaxation time   at last time step


# Initialize float and int constants to be of the correct type

if D == 3: beta_g /= 1.5                         # Allow for 3D

D           = np_int(D)                          # Number of dimensions
gamma_T     = np_float(gamma_T/dT_p)             # Exponent for temperature dependence
if D == 3: gamma_T /= 1.5                        # Allow in 3D for gamma_T applied to eta

NX          = np_int(nxi*nyi*nzi)                # Size of computational domain
_NX_        = np_int(nx_i*ny_i*nz_i)             # Size of communication domain

#  RKTLBM parameters 
#  Based on Reis and Phillips, 2007:
alfa_b    = 4/9
alfa_r    = 1.0 - (1.0 - alfa_b)
#  old method:
# alfa_r    = rho1/(rho1+rho2)
# alfa_b    = rho2/(rho1+rho2)

if not (0 <= alfa_r < 1):
   raise ValueError(f"Stability error: alfa_r={alfa_r} is out of bounds [0, 1).")

C_r    = np.zeros(na,dtype=np_float)
C_b    = np.zeros(na,dtype=np_float)

if D == 2:
    C_r[0] = alfa_r
    for a in range(1,5):  C_r[a] = (1.0-alfa_r)/5.0
    for a in range(5,na): C_r[a] = (1.0-alfa_r)/20.0
    C_b[0] = alfa_b
    for a in range(1,5):  C_b[a] = (1.0-alfa_b)/5.0
    for a in range(5,na): C_b[a] = (1.0-alfa_b)/20.0
elif D == 3:
    C_r[0] = alfa_r
    for a in range(1,7):  C_r[a] = (1.0-alfa_r)/12.0
    for a in range(7,na): C_r[a] = (1.0-alfa_r)/24.0
    C_b[0] = alfa_b
    for a in range(1,7):  C_b[a] = (1.0-alfa_b)/12.0
    for a in range(7,na): C_b[a] = (1.0-alfa_b)/24.0

rho1  = np_float(rho1)
rho2  = np_float(rho2)
rho12 = np_float(rho12)

# Convert Python strings to bytes which is equivalent to the C "char" type which
# can be passed to the C TLBM_fast() function

if   temperature == "reynolds":    Temperature = b"reynolds"
elif temperature == "arrhenius":   Temperature = b"arrhenius"
elif temperature == "unspecified": Temperature = b"unspecified"

if   nu_formula == "Tau":         Nu_formula = b"Tau"
elif nu_formula == "strain_rate": Nu_formula = b"strain_rate"

if   rheology == "linear":   Rheology = b"linear"
elif rheology == "archean":  Rheology = b"archean"
elif rheology == "archean2": Rheology = b"archean2"

# Normalize melt_position: strip whitespace and lowercase before matching.
# Otherwise a stray space in the .in file (e.g. "insitu ") leaves
# Melt_position undefined or unmatched in C → strcmp never returns 0
# → deposit code never runs and rho_r stays ≈ 0 even though ce_bs grows.
_melt_position_norm = str(melt_position).strip().lower()
if   _melt_position_norm == "surface": Melt_position = b"surface"
elif _melt_position_norm == "insitu" : Melt_position = b"insitu"
else:
    Melt_position = b""   # explicit empty-string sentinel; deposit will skip,
                          # but at least it won't carry stale stack memory.
    if rank == 0:
        print(f"WARNING: unrecognized melt_position={melt_position!r} "
              f"(normalized to {_melt_position_norm!r}). "
              f"Expected 'surface' or 'insitu'. Deposit will be skipped.",
              flush=True)
melt_position = _melt_position_norm   # keep Python side consistent too

if rank == 0: print(" %4.2f secs" %(time.time()-t0),flush=True)

# End   of array allocation
#==============================================================
# Initialize the get edges, BC and get_field
#==============================================================

if rank == 0: print("Initializing get_edges, BC, get_field & put_field ...",end="",flush=True)

t0 = time.time()

initialize_get_edges()

initialize_BC()

# Initialize fast get field

initialize_get_field_fast()
initialize_put_field_fast()

# Set an edge mask (used to avoid overflows in python code which calculates in edges)

if rank < nr:
    edge_mask = np.ones((ny_i,nz_i,nx_i),dtype=bool)
    edge_mask[yedge:yedge+my[Ry],zedge:zedge+mz[Rz],xedge:xedge+mx[Rx]] = False

if rank == 0: print(" %4.2f secs" %(time.time()-t0),flush=True)

#--------------------------------------------------------------
# End   initialization of get edges function on domains
#==============================================================
# Initialize variables
#==============================================================

if rank == 0: print("Initializing variables ...",flush=True)

#--------------------------------------------------------------
# Initialize T[z] and from this, the energy density and
# the thermal boundary condition arrays eta_u and eta_d
#--------------------------------------------------------------

if rank < nr:

    t0 = time.time()

    if rank == 0: print("    T and eta ...",end="",flush=True)

#   Initialize the thermal BC
    T_z = np.zeros(nz,dtype=np_float) # Temperature versus depth array
    P_z = np.zeros(nz,dtype=np_float) # Pressure versus depth array
    T_z[z_upper] = T_upper
    T_z[z_lower] = T_lower
    
    # Depth in km from the surface for each cell.
    # For annular geometry the cells are LOG-spaced in radius (not uniform),
    # so z_km must be computed from the actual radii array. For cartesian
    # geometry the cells are uniform and dx_p*z works.
    if geometry == "annular":
        # radii[0] = surface (R_p, largest), radii[nz-1] = CMB (smallest).
        # Depth from surface = R_p - radii[i] (in metres → km).
        z_km        = (R_p - radii) * 1e-3
        z_1d        = R_p - radii        # in metres (depth)
    else:
        z_1d        = dx_p * np.linspace(0,nz-1,nz)
        z_km        = z_1d * 1e-3
    if rank == 0:
        print("------------------------------------------------------------------",flush=True)
        print(f"Depth grid (km) — geometry={geometry}, nz={nz}, L_p={L_p*1e-3:.1f} km",flush=True)
        _Lp_km = L_p * 1e-3
        for _d in (0, 660, 2750, 2900):
            if _d > _Lp_km + 1e-3:        # skip targets deeper than the domain
                continue
            _i = int(np.argmin(np.abs(z_km - _d)))
            _label = "surface" if _d == 0 else ("CMB" if abs(_d - _Lp_km) < 1.0 else "")
            print(f"  z_km≈{_d:>4d} km  -> cell z={_i:>3d}  (actual {z_km[_i]:8.3f} km)  {_label}",flush=True)
        print("------------------------------------------------------------------",flush=True)

#   USERS should specify any temperature depth dependence below here
#   Namely, initialize T_z[z]

#   Safe default P0 (1 GPa everywhere) in case rheology doesn't set it
    P0 = put_field_z_2_xyz(np.ones(nz,dtype=np_float),default=1)

 
    if rheology == "archean":
       if T_z_file == "unspecified":
          # Branch A: BUILD all defaults (no files loaded).
          for z in range(nz):
#           Use thickness of lithosphere of 100 km or 120 km
            if z >= int(120 * 1000 * L_l/L_p):
                T_z[z] = T_lower * 0.70
            else:
                gradient_z = z/int(120 * 1000 * L_l/ L_p)
                T_z[z] = T_upper + gradient_z * (T_lower * 0.70 - T_upper)
          T_z[z_upper]  = T_upper
          T_z[z_lower]  = T_lower
          # T_z[z_lower] *= 1.5   # disabled 2026-11-16 (Amen): CMB = nominal T_lower = 3415 C
          # Pressure/Density for arrhenius equation
          rho_z    = rho_p * np.ones(nz,dtype=np_float)
          L_z      = np.zeros(nz,dtype=np_float)
          n_basalt = int(basalt_depth_p / dx_p)
          if n_basalt > 0:
              L_z[:n_basalt] = 1
          P_z  = np.zeros(nz)
          for i in range(nz):
              P_z[i] = z_km[i] * rho_z[i] * 9.81
          P_z *= 1e-6
          P0         = put_field_z_2_xyz(P_z,default=1) # Pressure in GPa

       elif T_z_file != "unspecified" and P_z_file == "unspecified":
          # Branch C: LOAD T_z only — build rho_z, L_z, P_z defaults.
          T_z_file   = run_dir + T_z_file + ".npy"
          T_z        = np_float(np.load(T_z_file))
          T_z       -= (T0_p + dT_p/2)
          rho_z      = rho_p * np.ones(nz, dtype=np_float)
          L_z        = np.zeros(nz, dtype=np_float)
          n_basalt   = int(basalt_depth_p / dx_p)
          if n_basalt > 0:
              L_z[:n_basalt] = 1
          P_z        = np.zeros(nz, dtype=np_float)
          for i in range(nz):
              P_z[i] = z_km[i] * rho_z[i] * 9.81
          P_z       *= 1e-6
          P0         = put_field_z_2_xyz(P_z, default=1) # Pressure in GPa

       else:
          # Branch D: LOAD everything (T_z and P_z files both specified).
          T_z_file   = run_dir + T_z_file   + ".npy"
          T_z        = np_float(np.load(T_z_file))
          T_z       -= (T0_p+dT_p/2)
          P_z_file   = run_dir + P_z_file   + ".npy"
          P_z        = np_float(np.load(P_z_file))
          P0         = put_field_z_2_xyz(P_z,default=1) # Pressure in GPa

    if RKTLBM :   
           RHO_r = np.zeros((ny,nz,nx),dtype=np_float)
           RHO_b = np.ones((ny,nz,nx),dtype=np_float)
           if rho_r_file == "fallingbox" or rho_b_file == "fallingbox":
              width_block  = nx//8
              height_block = nz//8
              top_block    = nz//8-height_block//8
              base_block   = top_block + height_block
              for y in range(ny):
                  for z in range(nz):
                      for x in range(nx):
                          RHO_r[y,z,x] = 0.0*rho1
                          RHO_b[y,z,x] = 1.0*rho2
                          if top_block<z<=base_block and nx//2-width_block//2<x<=nx//2+width_block//2:
                             RHO_r[y,z,x] = 1.0*rho1
                             RHO_b[y,z,x] = 0.0*rho2
           elif rho_r_file == "surfacebox" or rho_b_file == "surfacebox":
              # Wetting-angle test: put a TTG slab AT the upper wall.
              # Box spans middle half in x, depth 1..H in z (z=0 is the wall halo).
              width_block  = nx//4
              height_block = max(nz//16, 4)
              top_block    = 1                          # touch the upper wall (zu+1)
              base_block   = top_block + height_block
              for y in range(ny):
                  for z in range(nz):
                      for x in range(nx):
                          RHO_r[y,z,x] = 0.0*rho1
                          RHO_b[y,z,x] = 1.0*rho2
                          if (top_block <= z < base_block) and \
                             (nx//2 - width_block//2 <= x < nx//2 + width_block//2):
                              RHO_r[y,z,x] = 1.0*rho1
                              RHO_b[y,z,x] = 0.0*rho2
           elif rho_r_file == "unspecified" or rho_b_file == "unspecified":
                pass
           else:
               rho_r_file = run_dir + rho_r_file + ".npy"
               rho_b_file = run_dir + rho_b_file + ".npy"
               RHO_r      = np_float(np.load(rho_r_file))
               RHO_b      = np_float(np.load(rho_b_file))
 
           rho_r      = put_field(RHO_r)
           rho_b      = put_field(RHO_b)
           rh0_reff   = put_field(RHO_r+RHO_b)
           # Capture true t=0 mass sums BEFORE any kernel runs.
           if rank == 0:
               write_and_plot_field._buoy_init_r   = float(np.sum(RHO_r, dtype=np.float64))
               write_and_plot_field._buoy_init_b   = float(np.sum(RHO_b, dtype=np.float64))
               write_and_plot_field._buoy_init_tot = (write_and_plot_field._buoy_init_r
                                                    + write_and_plot_field._buoy_init_b)

    if rheology == "archean2":
        
        if T_z_file == "unspecified" and rho_z_file == "unspecified":
           # Default initial profiles: cold 120-km lithospheric lid, warm
           # adiabatic interior, hot bottom boundary layer.
           # Pure peridotite (L=0) everywhere — basalt evolves from melt.
           lith_z = max(int(120e3 * L_l / L_p), 1)
           for z in range(nz):
               if z >= lith_z:
                   T_z[z] = T_lower * 0.70
               else:
                   T_z[z] = T_upper + (z / lith_z) * (T_lower * 0.70 - T_upper)
           T_z[z_upper]  = T_upper
           T_z[z_lower]  = T_lower
           # T_z[z_lower] *= 1.5   # disabled 2026-11-16 (Amen): CMB = nominal T_lower = 3415 C

           rho_z    = rho2_p * np.ones(nz, dtype=np_float)
           L_z      = np.zeros(nz, dtype=np_float)
           n_basalt = int(basalt_depth_p / dx_p)
           if n_basalt > 0:
               L_z[:n_basalt] = 1                         # initial basalt surface layer
               rho_z[:n_basalt] = rho3_p 
           P_z   = np.zeros(nz, dtype=np_float)
           for i in range(nz):
               P_z[i] = z_km[i] * rho_z[i] * 9.81
           P_z *= 1e-6                                    # convert to GPa

           P0    = put_field_z_2_xyz(P_z,   default=1)    # Pressure in GPa
           rho_M = put_field_z_2_xyz(rho_z, default=1)
           L_ed  = put_field_z_2_xyz(L_z,   default=0)
           T_nd  = T_z
           T_nd  = put_field_z_2_xyz(np_float(T_nd))

        elif T_z_file == "test":
           for z in range(nz):
               if z >= int(120 * 1000 * L_l/L_p):
                  T_z[z] = T_lower * 0.70
               else:
                  gradient_z = z/int(120 * 1000 * L_l/ L_p)
                  T_z[z] = T_upper + gradient_z * (T_lower * 0.70 - T_upper)
           T_z[z_upper]  = T_upper
           T_z[z_lower]  = T_lower
           # T_z[z_lower] *= 1.5   # disabled 2026-11-16 (Amen): CMB = nominal T_lower = 3415 C

           rho_z = rho_p * np.ones(nz,dtype=np_float)
           L_ED   = np.zeros((ny,nz,nx),dtype=np_float)
           width_block  = nx//8
           height_block = nz//8
           top_block_r  = nz//8 - height_block//8
           base_block_r = top_block_r + height_block
           top_block_L  = nz - base_block_r
           base_block_L = nz - top_block_r
           for y in range(ny):
               for z in range(nz):
                   for x in range(nx):
                       L_ED[y,z,x] = 0.0
                       if top_block_L<z<=base_block_L and nx//2-width_block//2<x<=nx//2+width_block//2:
                          L_ED[y,z,x] = 1.0

           P_z  = np.zeros(nz)
           for i in range(nz):
               P_z[i] = z_km[i] * rho_z[i] * 9.81
           P_z *= 1e-6

           L_ed       = put_field(L_ED)
           P0         = put_field_z_2_xyz(P_z,default=1) # Pressure in GPa
           rho_M      = put_field_z_2_xyz(rho_z,default=1)
           T_nd       = T_z
           T_nd       = put_field_z_2_xyz(np_float(T_nd))

        elif T_z_file != "unspecified" and rho_z_file == "unspecified":
             T_z_file = run_dir + T_z_file   + ".npy"
             T_z      = np_float(np.load(T_z_file))
             rho_z    = rho2_p * np.ones(nz, dtype=np_float)
             L_z      = np.zeros(nz, dtype=np_float)
             n_basalt = int(basalt_depth_p / dx_p)
             if n_basalt > 0:
                L_z[:n_basalt] = 1                         # initial basalt surface layer
                rho_z[:n_basalt] = rho3_p  
             P_z   = np.zeros(nz, dtype=np_float)
             for i in range(nz):
                 P_z[i] = z_km[i] * rho_z[i] * 9.81
             P_z *= 1e-6                                    # convert to GPa

             P0    = put_field_z_2_xyz(P_z,   default=1)    # Pressure in GPa
             rho_M = put_field_z_2_xyz(rho_z, default=1)
             L_ed  = put_field_z_2_xyz(L_z,   default=0)
             T_nd  = (T_z - T0_p)/dT_p
             T_nd  = put_field_z_2_xyz(np_float(T_nd))
             T_z  -= (T0_p + dT_p/2)

        else:
           rho_z_file = run_dir + rho_z_file + ".npy"
           T_z_file   = run_dir + T_z_file   + ".npy"
           P_z_file   = run_dir + P_z_file   + ".npy"
           L_z_file   = run_dir + L_z_file   + ".npy"

           rho_z      = np_float(np.load(rho_z_file))
           P_z        = np_float(np.load(P_z_file))
           T_z        = np_float(np.load(T_z_file))      # Temperature in C
           L_z        = np_float(np.load(L_z_file))
           P0         = put_field_z_2_xyz(P_z,default=1) # Pressure in GPa
           rho_M      = put_field_z_2_xyz(rho_z,default=1)
           L_ed       = put_field_z_2_xyz(L_z,default=0)
           T_nd       = (T_z - T0_p)/dT_p
           T_nd       = put_field_z_2_xyz(np_float(T_nd))
           T_z       -= (T0_p+dT_p/2)
    z_nd            = z_1d/L_p
    z               = put_field_z_2_xyz(np_float(z_1d))
    z_nd            = put_field_z_2_xyz(np_float(z_nd))
    mask_pd         = np.where(z >  2e5 ,1,0).astype(bool)
    # Basalt -> TTG wherever the Moyen/Jain P-T windows permit (LP < 1 GPa,
    # MP 1-2.2 GPa, HP 2.35-5 GPa). The former 70 km floor was a band-aid for
    # the melt-factory runaway, which is now fixed properly (TVD depletion
    # transport, conservative budgets, exact population L source); it made the
    # run HP-only, contradicting the paper's Moyen 20/60/20 framing.
    mask_bs         = np.ones_like(z).astype(bool) 
    # Mantle Viscosity Profile (Arrhenius Equation) based on Jain et al. (2019)
    act_E           = np.zeros((nz),dtype=np_float)
    V_array         = np.zeros((nz),dtype=np_float)
    P_scale_array   = np.zeros((nz),dtype=np_float)
    Delta_eta_array = np.zeros((nz),dtype=np_float)
    for i in range(nz):
        if z_km[i] < 660:
           act_E[i]           = act_E_1
           V_array[i]         = act_V_1
           P_scale_array[i]   = P_scale_1
           Delta_eta_array[i] = Delta_eta_1
        elif z_km[i] < 2740:
           act_E[i]           = act_E_2
           V_array[i]         = act_V_2
           P_scale_array[i]   = P_scale_2
           Delta_eta_array[i] = Delta_eta_2
        else:
           act_E[i]           = act_E_3
           V_array[i]         = act_V_3
           P_scale_array[i]   = P_scale_3
           Delta_eta_array[i] = Delta_eta_3
    act_V           = (V_array * np.exp((-P_z*1e9)/P_scale_array)).astype(np_float)
    act_E           = put_field_z_2_z(act_E)
    act_V           = put_field_z_2_z(act_V)
    Delta_eta_1d    = Delta_eta_array.copy()          # keep (nz,) for nu_f_b depth profile
    Delta_eta_array = put_field_z_2_z(Delta_eta_array)

    # Set depth-dependent blue viscosity: only apply the 660/2740 km Delta_eta
    # profile for Arrhenius rheology.  For Reynolds the reference viscosity is
    # uniform (exp(-gamma_T*eta) handles the temperature dependence instead).
    if RKTLBM:
        if temperature == "arrhenius":
            nu_f_b = put_field_z_2_z(np_float(nu_f0_b * Delta_eta_1d), default=nu_f0_b)
        else:
            nu_f_b = put_field_z_2_z(np_float(nu_f0_b * Delta_eta_1d), default=nu_f0_b)

    if geometry == "annular":
       dx_rev = dx_ann[::-1] 
       dy_rev = dy_ann[::-1] 

#    if rank == 0:
#       print()
#       print("Input file check", flush=True)
#       for z in range(nz):
#           print(f"Temperature at depth {z} is T_z[{z}] = {T_z[z]:.2e}", flush=True)


#   Put the temperature into eta (energy density).
#   Build a global 3D temperature field T_xyz of shape (ny, nz, nx).
#   - If T_z came in as 1D depth profile  -> broadcast to (ny, nz, nx).
#   - If T_z came in as 3D from T_z_file  -> use directly.
#   Then distribute the 3D field to ranks via put_field (instead of the
#   1D-only put_field_z_2_xyz). This lets the user supply a fully 3D
#   T_z_file and have it placed cell-by-cell.
    if T_z.ndim == 1:
        T_xyz = np.broadcast_to(T_z[np.newaxis,:,np.newaxis], (ny,nz,nx)).copy()
    elif T_z.ndim == 3:
        T_xyz = np_float(T_z)
    else:
        if rank == 0:
            print(f"ERROR: T_z must be 1D (depth profile) or 3D (full field); got ndim={T_z.ndim}",
                  file=sys.stderr, flush=True)
        exit()
#   Cap the initial temperature at the peridotite solidus minus a margin so a
#   melt run does not start with a pre-molten deep layer. Wall rows keep their
#   BC values and a thin basal boundary layer remelts self-consistently.
    if rheology == "archean2" and melt_extraction:
        T_sol_init = np.where(P_z < 13.0, -5.1*P_z**2 + 133*P_z + 1120.85,
                              2045.0*(1.0 + P_z/92.0)**(1/1.3) - 273.15).astype(np_float)
        # The stored T field is offset from the melt module temperature by
        # T0_p + dT_p/2, so the cap must be expressed in the stored scale.
        _T_offset = T0_p + dT_p/2
        T_xyz[:,1:nz-1,:] = np.minimum(T_xyz[:,1:nz-1,:],
                                       (T_sol_init - 50.0 - _T_offset)[np.newaxis,1:nz-1,np.newaxis])
        # Conductive-lid cap: initial T limited to 10 K/km near the surface so
        # the initial crust starts below the 760 C TTG window floor (no t=0
        # conversion burst); deeper rows are unaffected (lid cap >> profile).
        _T_lid = 10.0 * z_km.astype(np_float)
        T_xyz[:,1:nz-1,:] = np.minimum(T_xyz[:,1:nz-1,:],
                                       (_T_lid - _T_offset)[np.newaxis,1:nz-1,np.newaxis])
        if rank == 0:
            print("Initial T capped at peridotite solidus - 50 K below the surface rows",flush=True)

    eta  = put_field(T_xyz, default=0)
    del T_xyz

#   Add sinusoudal temperature fluctuations

#   if amp_sinusoid > 0 and rank < nr: eta = initialize_sinusoidal_temperature(eta)

#   Add random temperature fluctuations throughput the model

    eta[:,:,:] += np.random.normal(0.0,T_p_random,size=(ny_i,nz_i,nx_i))

#   Add random temperature fluctuations on the boundaries to initiate convection

    np.random.seed(rank)
    if boundary_u:
        eta[:,zu,:] += np.random.normal(0.0,dT_p_random,size=(ny_i,nx_i))
    if boundary_d:
        eta[:,zd,:] += np.random.normal(0.0,dT_p_random,size=(ny_i,nx_i))

    del T_z

#   Add a sinusoidal fluctuation to initiate convection if initialize_ncells > 0

    if initialize_ncells > 0 and rank < nr:
        if geometry == "annular" or (geometry == "spherical" and D == 2):
            dtheta             = np.pi*(theta1/180+1/initialize_ncells)
            initialize_ncells *= (theta2-theta1)/360
        else:
            dtheta             = 0
        x0 = Rx * nxi
        x1 = x0 + mx[Rx]
        dT_sinusoid = np.zeros(nx_i,dtype=np_float)
        x_sinusoid  = np.linspace(x0,x1,nxi,dtype=np_float)
        sinusoidal_fluct = np.cos(2*np.pi*x_sinusoid*initialize_ncells/nx+dtheta)
        power_sinusoid   = 10
        sinusoidal_fluct = np.sign(sinusoidal_fluct)*sinusoidal_fluct**power_sinusoid
        dT_sinusoid[xedge:xedge+nxi] = cell_amplitude*np.cos(2*np.pi*x_sinusoid*initialize_ncells/nx+dtheta)
        if boundary_u: eta[:,zu,xedge:xedge+nxi] += dT_sinusoid[xedge:xedge+nxi]
        if boundary_d: eta[:,zd,xedge:xedge+nxi] += dT_sinusoid[xedge:xedge+nxi]
        del dT_sinusoid

#   Set eta_u and eta_d which are used to impose the
#   thermal BC on the upper and lower boundaries

#   Optional 2D lithology IC from file (mirrors rho_r_file): L_file = <name>
#   loads run_dir/<name>.npy with shape (ny,nz,nx), values in [0,1]; archean2 only.
    if L_file != "unspecified":
        if rheology == "archean2":
            _L_path = run_dir + L_file + ".npy"
            L_ed    = put_field(np_float(np.load(_L_path)))
            if rank == 0: print("Loaded 2D lithology IC from %s" % (_L_path), flush=True)
        elif rank == 0:
            print("L_file specified but rheology != archean2: IGNORED", flush=True)

#   Optional 2D lithology IC from file (mirrors rho_r_file): L_file = <name>
#   loads run_dir/<name>.npy with shape (ny,nz,nx), values in [0,1]; archean2 only.
    if L_file != "unspecified":
        if rheology == "archean2":
            _L_path = run_dir + L_file + ".npy"
            L_ed    = put_field(np_float(np.load(_L_path)))
            if rank == 0: print("Loaded 2D lithology IC from %s" % (_L_path), flush=True)
        elif rank == 0:
            print("L_file specified but rheology != archean2: IGNORED", flush=True)

    if boundary_u: eta_u[:,:] = eta[:,zu,:]
    if boundary_d: eta_d[:,:] = eta[:,zd,:]

# Convert temperature to energy density

    eta *= DR_2

    if rank == 0: print(" %4.2f secs" %(time.time()-t0),flush=True)
#--------------------------------------------------------------
#   Initialize nu_f and nu_g which can be depth dependent
#--------------------------------------------------------------

    t0 = time.time()
    if rank == 0: print("    nu_f and nu_g ...",end="",flush=True)
    
#   These arrays allow nu_f and nu_g to be depth dependent

    NU_F   = nu_f0*np.ones(nz,dtype=np_float)
    NU_G   = nu_g0*np.ones(nz,dtype=np_float)
    nu_f0_initial = nu_f0

#   USERS should specify any viscosity or thermal diffusivity depth dependence below here
#   Namely, initialize NU_F[z] and NU_G[z]

#   USERS should specify any viscosity or thermal diffusivity depth dependence above here

#   Rescale NU_F and NU_G to allow for spherical geometry if needed

    if ny == 1 and geometry == "spherical": geometry = "annular"

    if geometry == "annular":
        NU_F   *= nu_factor
        NU_G   *= nu_factor
        if RKTLBM:
            if temperature == "arrhenius":
                nu_f_b = put_field_z_2_z(np_float(nu_f0_b * nu_factor * Delta_eta_1d), default=nu_f0_b)
            else:
                nu_f_b = put_field_z_2_z(np_float(nu_f0_b * nu_factor * Delta_eta_1d), default=nu_f0_b)
            nu_f_r = put_field_z_2_z(np_float(nu_f0_r * nu_factor), default=nu_f0_r)

    if code == "python":
        nu_f   = put_field_z_2_xyz(NU_F,default=nu_f0)
        nu_g   = put_field_z_2_xyz(NU_G,default=nu_g0)
        
    elif code == "c" and not RKTLBM:
        nu_f   = put_field_z_2_z(NU_F,default=nu_f0)
        nu_g   = put_field_z_2_z(NU_G,default=nu_g0)

    if not RKTLBM:
       tau_f   = nu_f * 3 + 0.5
       tau_g   = nu_g * 3 + 0.5
       tau_h   = nu_g * Pr/Peclet * 3 + 0.5
    else:
       tau_h   = nu_f_b / Peclet * 3 + 0.5

    if code == "python":
        nu_f_last [:,:,:] = nu_f [:,:,:]
        tau_f_last[:,:,:] = tau_f[:,:,:]
    elif code == "c" and not RKTLBM:
        for z in range(nz_i):
            nu_f_last [:,z,:] = nu_f [z]
            tau_f_last[:,z,:] = tau_f[z]
    elif code == "c" and RKTLBM:
        for z in range(nz_i):
            nu_f_last [:,z,:] = nu_f_b[z]
            tau_f_last[:,z,:] = nu_f_b[z] * 3 + 0.5

    del NU_F, NU_G

#   g ~ dx^3 annular gauge (benchmark test_option S): per-row Boussinesq factor
    BETA_G_NZ = beta_g * np.ones(nz, dtype=np_float)
    if g_dx_gauge and geometry == "annular":
        BETA_G_NZ = np_float(BETA_G_NZ * (radii/radii[0])**3)
        if rank == 0:
            print("g_dx_gauge ACTIVE: beta_g scaled by (dx(z)/dx(0))^3, surface %.4g -> base %.4g" % (BETA_G_NZ[0], BETA_G_NZ[-1]), flush=True)
    elif g_dx_gauge:
        if rank == 0: print("g_dx_gauge requested but geometry != annular: IGNORED", flush=True)
    if code == "c":
        BETA_G = put_field_z_2_z(BETA_G_NZ, default=beta_g)
    else:
        BETA_G = np_float(BETA_G_NZ + 0)


    if rank == 0: print(" %4.2f secs" %(time.time()-t0),flush=True)

#--------------------------------------------------------------
#   Initialize the plastic rheology which can be depth dependent
#--------------------------------------------------------------
    if rank == 0: print("    Tau_y, N and Delta_g ...",end="",flush=True)

#   Set default values

    if rheology == "linear":
        Tau_y_default = 1
        N_default     = 0
    else:
        Tau_y_default = 1.e30
        N_default     = N

#   Declare the array for radiogenic heat generation

    Delta_g = np.zeros(na,dtype=np_float)

#   Declare arrays for Tau y and N as a function of z

    TAU_Y   = Tau_y_default*np.ones(nz,dtype=np_float)   # Yield stress
    NN      = N_default    *np.ones(nz,dtype=np_float)   # Power exponent

#   USERS should specify any yield stress or exponent depth dependence below here
#   Namely, initialize TAU_Y[z] and NN[z]

#   Archean yield stress versus depth

    if rheology == "archean" or rheology == "archean2":
       if rank == 0:
          print() 
          print("TAU physical from input file")
       for z in range(nz):
           TAU_Y[z] = (Tau_y0_p + (z * dx_p) * Tau_yi_p)
            # print(f"tau physical at depth {z} is TAU_Y[z] = {Tau_y0_p + (z * dx_p) * Tau_yi_p:.2e}", flush=True) 

#       Convert TAU_Y from physical to lattice units which uses the following:
#           - density in kg/m^3 to mass per lattice-site (1/rho_p)
#           - acceleration in m/s^2 to lattice units     (g_l/g_p)
#           - force to stress in lattice units           (1/dx_p**(D-1))

       if   rheology == "archean" : conversion_Tau_p_2_l = (1/rho_p)
       elif rheology == "archean2": conversion_Tau_p_2_l = (1/rho_z)
 
       conversion_Tau_p_2_l *= (g_l/g_p) / dx_p**(D-1)
       TAU_Y                *= conversion_Tau_p_2_l

#        if rank == 0:
#           print("TAU lattice calculated")
#           for z in range(nz):
#               print(f"tau lattice at depth {z} is TAU_Y[{z}] = {TAU_Y[z]:.2e}", flush=True)       

#   USERS should specify any yield stress or exponent depth dependence above here

#   Put the specified depth dependence of the yield stress and
#   exponent into the Tau_y and N arrays needed by python and c

    if code == "python":

        Tau_y         = put_field_z_2_xyz(TAU_Y,default=Tau_y_default)
        N             = put_field_z_2_xyz(NN   ,default=N_default    )
        nu_star       = np_float(nu_star) 
    elif code == "c":
        Tau_y         = put_field_z_2_z  (TAU_Y,default=Tau_y_default)
        N             = put_field_z_2_z  (NN   ,default=N_default    )
        nu_star       = np_float(nu_star)        

    del NN, TAU_Y

    if rank == 0: print(" %4.2f secs" %(time.time()-t0),flush=True)

#==============================================================
#   Initialize distributions
#==============================================================

t0 = time.time()

if rank == 0: print("Initializing distributions ...",end="",flush=True)

if rank < nr:

#   Initialize distributions

    initialize_distributions_fast()     # Initialize f and g distributions

    # Reference density arrays for RKTLBM shifted distribution.
    # The reference MUST be spatially uniform (a constant scalar broadcast to 3-D)
    # so that rho_r_init[source] = rho_r_init[dest] during streaming and the MAC
    # error term (rho_r_init[dest] − Σ_a rho_r_init[src_a]·w[a]) is identically 0.
    # Using rho_r.copy() (the spatially-varying initial field) is WRONG: cells at
    # the moving interface would have rho_r_init[source] ≠ rho_r_init[dest], giving
    # a persistent MAC error proportional to the density jump that grows as the
    # interface sweeps through the domain.
    # Always allocate when RKTLBM=True so ctypes always has a valid pointer to pass;
    # the C kernel only reads these arrays when shifted_distribution=True.
    if RKTLBM:
        if shifted_distribution:
            rho_r_init = np.full(rho_r.shape, rho1, dtype=np_float)
            rho_b_init = np.full(rho_b.shape, rho2, dtype=np_float)
        else:
            rho_r_init = np.zeros(rho_r.shape, dtype=np_float)   # ignored by C
            rho_b_init = np.zeros(rho_b.shape, dtype=np_float)   # ignored by C
    else:
        rho_r_init = None
        rho_b_init = None

    if rank == 0: print(" %4.2f secs" %(time.time()-t0),flush=True)

#==============================================================
# Initialize streaming
#==============================================================

#--------------------------------------------------------------------------
# Initialize streaming index arrays
#--------------------------------------------------------------------------

t0 = time.time()

if rank == 0: print("Initializing streaming index arrays ...",end="",flush=True)

#   Initialize arrays for x-x_a, y-y_a and z-z_a

xx   = np.linspace(0,nx_i-1,nx_i,dtype=np_int)
yy   = np.linspace(0,ny_i-1,ny_i,dtype=np_int)
zz   = np.linspace(0,nz_i-1,nz_i,dtype=np_int)
x_xa = np.zeros((na,nx_i),dtype=np_int)
y_ya = np.zeros((na,ny_i),dtype=np_int)
z_za = np.zeros((na,nz_i),dtype=np_int)

for a in range(na):
    x_xa[a,xx] = (xx-c[a,X]+nx_i)%nx_i
    y_ya[a,yy] = (yy-c[a,Y]+ny_i)%ny_i
    z_za[a,zz] = (zz-c[a,Z]+nz_i)%nz_i

#--------------------------------------------------------------------------
# Initialize vector streaming arrays for python if needed
#--------------------------------------------------------------------------

if code == "python":

#   Indices for vectorized streaming

    _x_xa_ = np.zeros((na,ny_i*nz_i*nx_i),dtype=np_int)

    xx_xa  = np.zeros((na,ny_i,nz_i,nx_i),dtype=np_int)
    yy_ya  = np.zeros((na,ny_i,nz_i,nx_i),dtype=np_int)
    zz_za  = np.zeros((na,ny_i,nz_i,nx_i),dtype=np_int)

#   Note: meshgrid has a strange default indexing = "xy" which transposes
#   the first two axes. Use of indexing = "ij" avoids this problem

    for a in range(na):
        yya,zza,xxa = np.meshgrid(y_ya[a],z_za[a],x_xa[a],indexing="ij")
        xx_xa[a] = xxa
        yy_ya[a] = yya
        zz_za[a] = zza
    del xxa, yya, zza

    xx_xa  = xx_xa.reshape(na,_NX_) ; yy_ya = yy_ya.reshape(na,_NX_) ; zz_za = zz_za.reshape(na,_NX_)
    for a in range(na):
        _x_xa_[a]     = yy_ya[a]*nz_i*nx_i + zz_za[a]*nx_i + xx_xa[a]

    del xx_xa, yy_ya, zz_za

if rank == 0: print(" %4.2f secs" %(time.time()-t0),flush=True)

#==========================================================================
# Initialize variables for TLBM calculations time loop
#==========================================================================
#   Flops per lattice site per time step calculation
#--------------------------------------------------------------------------

Flops = 0

if shifted_distribution: Flops += 1 # rho = 1 + drho
if code == "c":
    Flops += (na-1)*2+1+(na-1)*D+D
else:
    Flops += na*2+1+na*D+D
Flops += 1 + 2*D-1               # u2
Flops += 1 + (na-1)*(2*D-1+1+5)  # cu and equil core (1+3cu+(9/2)cu^2-(3/2)u2)
if shifted_distribution:
    Flops += 2                   # Allow for + drho and + w
Flops += na*9                    # Delta_f/g + collision
if D == 2:                       # Boussinesq
    Flops += 2 + 6*2
else:
    Flops += 2 + 10*2
if gamma_T > 0:
    Flops += 2
if rheology != "linear":
    Flops += (D**2*(na-1)*2+2+D**2+2*D**2+1+2)
    Flops += (4+4)

#--------------------------------------------------------------------------
# Initialize timing variables for various parts of the algorithm
#--------------------------------------------------------------------------

totalTimeEdges      = zero; # Total time to get edges
totalTimeEdges1     = zero; # Time to initiate send/recv requests               (overlap case)
totalTimeEdges2     = zero; # Time to finalize send/recv requests & get corners (overlap case)
totalTimeStream     = zero; # Time of streaming
totalTimeCalcs      = zero; # Total time for entire LBM calculations
totalTimeCalcs1     = zero; # Time of calculations (interior regions - overlap case)
totalTimeCalcs2     = zero; # Time of calculations (edge     regions - overlap case) 
totalTimeMach       = zero; # Time of mach number check
totalTimeOutput     = zero; # Time of output

time0               = time.time() # Set initial time for elapsed time calculation

#--------------------------------------------------------------------------
# Initialize the computed time step t and integral dt = LBM time
# and int(integral dt) and the last output/print/timing which is
# used to ensure that there is a single output/print/timing even
# for dt < 1 when several values of t have int(integral dt)%Inc=0
 
# Initialize integer LBM time for outputs

it_last_data       = -100000000000
it_last_plots      = -100000000000
it_last_timing     = -100000000000
it_last_mach       = -100000000000
it_last_vary_dt    = -100000000000
 
t_computed =  0       # Initialize computed time steps
t          =  0       # Initialize LBM time = integral dt where dt may be variable
dt         =  1       # Initial    LBM timestep
dt_last    =  2       # Initial    LBM timestep

mach_overflow = False                # Used to exit if mach_max is NAN, INF, or > mach_limit_to_exit
dT_overflow   = False                # Used to exit if mach_max is NAN, INF, or > dT_limit_to_exit
dT_max        = dT_p+2*(dT_p_random+cell_amplitude)
dT_limit_to_exit = 10*dT_max

factor        = 1                    # Initial factor for timestep
factor_max    = factor_max_vary_dt   # Max timestep factor
# Adaptive dt floor. Cutting dt rescales nu_g toward zero numerical
# diffusivity and the thermal AD field destabilizes once nu_g*factor drops
# below ~1e-5 lattice (M1_hd_smC blew up at nu_g*factor ~ 8e-6, held at
# 1.7e-5). A hard floor of 1.0 is safe but forbids the controller from
# arresting fast basal jets (M1_hpc: mach crept to 0.21 at the floor and the
# lattice nan-ed at 110 My). Allow cuts down to the thermal margin, never
# below 0.2 of base.
factor_min    = max(0.1, 2.5e-5/max(float(nu_g0), 1.0e-30))  # 2x margin: smH blew at nu_g*dt=1.1e-5
factor_min    = min(factor_min, 1.0) # never ABOVE base either
if rank == 0:
    print("vary-dt floor: factor_min = %6.3f (thermal-margin adaptive)" % factor_min, flush=True)

# Do not vary dt if factor_max_vary_dt is 1 or less
if abs(factor_max_vary_dt - 1) < 1.e-3: start_vary_dt = nt

it_mach_calc  = -1                   # Last mach calulation
it_output     = -1                   # Last output time
rho_bar       =  1                   # Initial lattice density mean
new_rho_bar   = False                # Set new rho_bar  bool
new_mach_max  = False                # Set new mach_max bool
mach_max_run    = 0                  # Max              mach for run
it_mach_max_run = 0                  # Time step of max mach for run
mach_max_sum    = 0                  # Sum       of   max mach for run
n_mach_max_sum  = 0                  # Number    of   max mach summed (to compute average)
mach_max_ave    = 0                  # Average   of   max mach for run

#=======================================================================================
# Start of time loop to do TLBM calculations
#=======================================================================================

if rank == 0:
    print("==================================================================",flush=True)
    now = datetime.now()
    date_string =  now.strftime("%B %d, %Y, %H:%M:%S")
    print("Starting TLBM simulation: %s" %(date_string),flush=True)
    print("==================================================================",flush=True)

# --- MASS-SOURCE DIAGNOSTIC: print wetting BC fractions so we know
#     if halos are injecting red fluid every step ---
# if rank == 0 and RKTLBM:
#     try:
#         print(f"[WETTING BCs] BC_r_wetting_u={BC_r_wetting_u:.6f}  BC_b_wetting_u={BC_b_wetting_u:.6f}  "
#               f"BC_r_wetting_d={BC_r_wetting_d:.6f}  BC_b_wetting_d={BC_b_wetting_d:.6f}  "
#               f"boundary_u={boundary_u}  boundary_d={boundary_d}  "
#               f"boundary_l={boundary_l}  boundary_r={boundary_r}",
#               flush=True)
#         print(f"[DOMAIN] nz_i={nz_i} zedge={zedge} zu={zu} zd={zd} "
#               f"ny_i={ny_i} yedge={yedge} nx_i={nx_i} xedge={xedge}",
#               flush=True)
#     except Exception:
#         pass

# Surface melt deposition is handled entirely by the C kernels:
#   * RKTLBM_rest deposits melt into rho_r[zu] and recolors f_r[zu]
#     so that sum(f_r[zu][a]) == rho_r[zu] (see TLBM_rest_kernel.c).
#   * RKTLBM_strmac streams f_r conservatively and recomputes rho_r
#     from the post-streaming f_r sum.
# No Python save/restore is needed; doing so would inject phantom
# mass each timestep (saved value re-imposed on top of mass that
# already streamed into interior cells -> exponential growth of rho_r).


# def _dbg_redsum(tag):
#     _R = get_field(rho_r)
#     if rank == 0:
#         print("DBGRED %-6s it=%d sum_r=%.9e" % (tag, it, float(np.sum(_R, dtype=np.float64))), flush=True)

while int(t-dt) <= nt:

#   With variable dt, t and dt must be bitwise identical on every rank: it is
#   derived from int(t*dt_dt0), a knife-edge truncation, and a sub-ULP per-rank
#   difference makes ranks take different Inc_* branches, which deadlocks the
#   collectives inside them (M1_hd_smC hung twice this way). Broadcast the two
#   scalars from rank 0 once per step; the cost is microseconds.
    _tdt = comm.bcast((t, dt), root=0)
    t, dt = _tdt

#   Global TTG budget: while the cumulative TTG DI (rank 0's latest growth
#   diagnostic) is at or above TTG_DI_max, basalt->TTG extraction is paused by
#   making the melt-window trigger M2_bs unreachable. Re-enables automatically
#   if DI falls back below the cap.
    if rank == 0:
        _cap_state = (len(DI_bs) > 0 and DI_bs[-1] >= TTG_DI_max)
        _tot = DI_bs[-1] if len(DI_bs) > 0 else 0.0
        if _tot > 0.05 and len(DI_bs_LP) > 0:
            _gates = []
            for _sh, _tgt in ((DI_bs_LP[-1]/_tot, TTG_LP_frac),
                              (DI_bs_MP[-1]/_tot, TTG_MP_frac),
                              (DI_bs_HP[-1]/_tot, TTG_HP_frac)):
                _gates.append(1.0 if _sh <= _tgt else 0.0)   # hard quota: threshold-driven extraction defeats fractional gates
        else:
            _gates = [1.0, 1.0, 1.0]
        _caps = (_cap_state, _gates[0], _gates[1], _gates[2])
    else:
        _caps = None
    _caps = comm.bcast(_caps, root=0)
    if rank == 0:
        if _caps[0] and not bs_capped:
            print("TTG cap engaged: DI_bs = %.4f %% >= TTG_DI_max = %.4f %% - basalt extraction paused" % (DI_bs[-1], TTG_DI_max), flush=True)
        for _wn, _g, _was in (("LP", _caps[1], lp_gatev), ("MP", _caps[2], mp_gatev), ("HP", _caps[3], hp_gatev)):
            if _g < 0.999 and _was >= 0.999:
                print("TTG %s window throttled: gate = %.3f (share above target, proportional)" % (_wn, _g), flush=True)
            if _g >= 0.999 and _was < 0.999:
                print("TTG %s window throttle released" % _wn, flush=True)
    bs_capped, lp_gatev, mp_gatev, hp_gatev = _caps

    it = int(t*dt_dt0) # Calculate the integer time. Note that the division by dt_dt0
    # _dbg_on = comm.bcast(bool((2000 <= it <= 2020) or (2245 <= it <= 2262) or (30000 <= it <= 30030) or (60000 <= it <= 60030) or (100000 <= it <= 100030) or ((it % 2000) < 3)), root=0)
    # if _dbg_on: _dbg_redsum("preK")
                       # ensures the integer time is relative to the specified physical
                       # time step, even if nu_l had to be adjusted for the nu_l_min constraint

#   Get the edges of domains

    if nr > 1 and rank < nr:
        t0 = time.time()
        if code == "c" and overlap and not RKTLBM:
            get_edges_start() # Start send/recv requests for edges
        else:
            get_edges()       # Get all edges
        if overlap: totalTimeEdges1 += (time.time()-t0)
        else:       totalTimeEdges  += (time.time()-t0)

#   Do TLBM calculations only on active ranks

    if rank < nr:

#       Swap pointers of f, g and h (ie. one-step streaming) to minimize the cost of streaming
        if not RKTLBM:
           f_tmp = f ; f = f_last ; f_last = f_tmp
           g_tmp = g ; g = g_last ; g_last = g_tmp
           h_tmp = h ; h = h_last ; h_last = h_tmp      

#       USERS can change as required the function heat_from_radiogenic_decay()

#       Allow for heat from radiogenic decay

        if rheology == "archean" or rheology == "archean2":
            Delta_g[:] = heat_from_radiogenic_decay(t)*w[:]

#-------------------------------------------------
#       C code option
#-------------------------------------------------
        if code == "c" and RKTLBM:
                t0 = time.time()

                # hd/hb now ride the main packed edge exchange (slots i_HD/i_HB,
                # face-crossing direction subsets, persistent requests) instead of
                # the standalone blocking get_edges_hd() exchange.

                RKTLBM_strmac( \
                na,D,nx_i,ny_i,nz_i,xedge,yedge,zedge, \
                mx_i,my_i,mz_i, \
                g, g_last, \
                x_xa,y_ya,z_za, \
                rho_r, rho_b, \
                f_r, f_r_last, \
                f_b, f_b_last, \
                Rheology, \
                h, h_last, \
                hd, hd_last, \
                hb, hb_last, \
                boundary_l,xl,naL,nbL,bL,biL, \
                boundary_r,xr,naR,nbR,bR,biR, \
                boundary_b,yb,naB,nbB,bB,biB, \
                boundary_f,yf,naF,nbF,bF,biF, \
                boundary_u,zu,naU,nbU,bU,biU, \
                boundary_d,zd,naD,nbD,bD,biD, \
                BC_r_wetting_u, BC_b_wetting_u, BC_r_wetting_d, BC_b_wetting_d, \
                shifted_distribution, \
                rho_r_init, rho_b_init, \
                rank)

                # Melt extraction at certain time in Myr True or False
                if np_int(it*dt_p/sec_my) >= melt_time :
                   melt_extraction = True
                else:
                   melt_extraction = False
                # Melt extraction position surface, insitu or False
                # Surface deposit is now handled in C (TLBM_rest_kernel.c, in the
                # A2_MELT block) so the python deposit below is a no-op. The C
                # version updates rho_r/rho_b AND the underlying f_r/f_b
                # distributions at z=zu, so the deposit survives the next strmac
                # rho_r=sum(f_r) recomputation. Re-enable the python block below
                # if you ever remove the C-side surface deposit.
                if rheology == "archean2" and melt_extraction:
                   if melt_position == "surface":
                      pass
                      # --- old python surface deposit (now handled in C) ---
                      # M_BS         = get_field(M_ext_bs)
                      # RHO_b        = get_field(rho_b)
                      # RHO_r        = get_field(rho_r)
                      # DEPOSIT      = np.zeros((ny,nz,nx),dtype=np_float)
                      # if rank == 0:
                      #    DEPOSIT[:,0,:]  = np.sum(M_BS, axis=1)
                      #    RHO_b          -= DEPOSIT
                      #    RHO_r          += DEPOSIT
                      # rho_b = put_field(RHO_b)
                      # rho_r = put_field(RHO_r)

                   else:
                      pass

                if nr > 1 and rank < nr:
                   get_edges_rho()

                RKTLBM_rest( \
                na,D,nx_i,ny_i,nz_i,xedge,yedge,zedge, \
                mx_i,my_i,mz_i, \
                w, c, \
                eta, \
                u, \
                g_last, \
                x_xa,y_ya,z_za, \
                beta_g, \
                tau_f0_r, tau_f0_b, tau_g, \
                rho_r, rho_b, \
                f_r_last, \
                f_b_last, \
                C_r, C_b, \
                rho1_p, rho2_p, rho3_p, beta_g_factor, g_comp_factor, \
                FF, \
                A, B, interface_thickness, g_comp, \
                BC_r_wetting_u, BC_b_wetting_u, BC_r_wetting_d, BC_b_wetting_d, \
                Rheology, Temperature, Nu_formula,\
                tau_h, \
                nu_f_r, nu_f_b, nu_star, \
                nu_l_min, nu_l_min_cap, nu_l_max_cap, gamma_T, Tau_y, N, Delta_g, \
                nu_total, tau_f_last, \
                T_nd, \
                act_E, act_V, Delta_eta_array, eta_p_input, DR_2, T0_p, dT_p, z_nd, \
                h_last, L_ed, \
                hd_last, D_ed, \
                hb_last, D_eb, \
                P0, cumulative_extracted_pd, cumulative_extracted_bs, \
                cumulative_extracted_LP, cumulative_extracted_MP, cumulative_extracted_HP, \
                P1_thresh, P2_thresh, \
                M1_pd, M2_pd, M4_pd, \
                M1_bs, (np_float(2.0) if bs_capped else M2_bs), M4_bs, \
                collision, mask_bs, mask_pd, \
                M_ext_bs_LP_TTG, M_ext_bs_MP_TTG, M_ext_bs_HP_TTG,
                rho0_pd, rho_melt_pd, M_ext_pd, \
                rho0_bs, rho_melt_bs, \
                rho_solid, rho_eff, \
                g_l, rho_M, rho_reff, \
                melt_extraction, Melt_position, \
                boundary_l,xl,naL,nbL,bL,biL, \
                boundary_r,xr,naR,nbR,bR,biR, \
                boundary_b,yb,naB,nbB,bB,biB, \
                boundary_f,yf,naF,nbF,bF,biF, \
                boundary_u,zu,naU,nbU,bU,biU, \
                boundary_d,zd,naD,nbD,bD,biD, \
                eta_u,eta_d, \
                shifted_distribution, \
                rho_r_init, rho_b_init, \
                rank, \
                np_float(lp_gatev), np_float(mp_gatev), np_float(hp_gatev))

                f_r_tmp = f_r ; f_r = f_r_last ; f_r_last = f_r_tmp
                f_b_tmp = f_b ; f_b = f_b_last ; f_b_last = f_b_tmp
                g_tmp   = g   ; g   = g_last   ; g_last   = g_tmp
                h_tmp   = h   ; h   = h_last   ; h_last   = h_tmp
                hd_tmp  = hd  ; hd  = hd_last  ; hd_last  = hd_tmp
                hb_tmp  = hb  ; hb  = hb_last  ; hb_last  = hb_tmp

                totalTimeCalcs += (time.time()-t0)
                # if _dbg_on: _dbg_redsum("postK")

        if code == "c" and not RKTLBM:

#-------------------------------------------------
#           C calculations, case of no overlapping
#-------------------------------------------------

            if not overlap or nr == 1:
            
                t0 = time.time()

#               Calculations in all of domain

                overlap = False ; interior = True

                TLBM_fast(na,D,nx_i,ny_i,nz_i,mx_i,my_i,mz_i,xedge,yedge,zedge, \
                w,c,shifted_distribution,overlap,interior, \
                rho,eta,u, \
                f_last,f,g_last,g, \
                x_xa,y_ya,z_za, \
                nu_f, tau_g, tau_h, BETA_G, \
                Temperature, Rheology, Nu_formula, \
                nu_l_min, nu_l_min_cap, nu_l_max_cap, gamma_T, Tau_y, N, Delta_g, \
                nu_f_last, tau_f_last, nu_total, \
                nu_star, \
                P1_thresh, P2_thresh, \
                M1_pd, M2_pd, M3_pd, M4_pd, \
                M1_bs, (np_float(2.0) if bs_capped else M2_bs), M3_bs, M4_bs, \
                rho0_pd, rho_melt_pd, M_ext_pd, \
                rho0_bs, rho_melt_bs, M_ext_bs, \
                M_ext_bs_LP_TTG, M_ext_bs_MP_TTG, M_ext_bs_HP_TTG, \
                mask_bs, mask_pd, \
                g_l, \
                T_nd, \
                cumulative_extracted_pd, cumulative_extracted_bs, \
                cumulative_extracted_LP, cumulative_extracted_MP, cumulative_extracted_HP, \
                P0, rho_M, L_ed,\
                h_last, h, \
                DR_2, T0_p, dT_p, z_nd, \
                rho_solid, rho_eff, \
                boundary_l,xl,naL,nbL,bL,biL,boundary_r,xr,naR,nbR,bR,biR, \
                boundary_b,yb,naB,nbB,bB,biB,boundary_f,yf,naF,nbF,bF,biF, \
                boundary_u,zu,naU,nbU,bU,biU,boundary_d,zd,naD,nbD,bD,biD, \
                eta_u,eta_d, \
                pre_streaming_BC, \
                rank)

                totalTimeCalcs += (time.time()-t0)

#-------------------------------------------------
#           C calculations, case of overlapping
#-------------------------------------------------

            else:

#               Calculations in interior regions of domains

                t0 = time.time()

                overlap  = True ; interior = True

                TLBM_fast(na,D,nx_i,ny_i,nz_i,mx_i,my_i,mz_i,xedge,yedge,zedge, \
                w,c,shifted_distribution,overlap,interior, \
                rho,eta,u, \
                f_last,f,g_last,g, \
                x_xa,y_ya,z_za, \
                nu_f, tau_g, tau_h, BETA_G, \
                Temperature, Rheology, Nu_formula, \
                nu_l_min, nu_l_min_cap, nu_l_max_cap, gamma_T, Tau_y, N, Delta_g, \
                nu_f_last, tau_f_last, nu_total, \
                nu_star, \
                P1_thresh, P2_thresh, \
                M1_pd, M2_pd, M3_pd, M4_pd, \
                M1_bs, (np_float(2.0) if bs_capped else M2_bs), M3_bs, M4_bs, \
                rho0_pd, rho_melt_pd, M_ext_pd, \
                rho0_bs, rho_melt_bs, M_ext_bs, \
                M_ext_bs_LP_TTG, M_ext_bs_MP_TTG, M_ext_bs_HP_TTG, \
                mask_bs, mask_pd, \
                g_l, \
                T_nd, \
                cumulative_extracted_pd, cumulative_extracted_bs, \
                cumulative_extracted_LP, cumulative_extracted_MP, cumulative_extracted_HP, \
                P0, rho_M, L_ed,\
                h_last, h, \
                DR_2, T0_p, dT_p, z_nd, \
                rho_solid, rho_eff, \
                boundary_l,xl,naL,nbL,bL,biL,boundary_r,xr,naR,nbR,bR,biR, \
                boundary_b,yb,naB,nbB,bB,biB,boundary_f,yf,naF,nbF,bF,biF, \
                boundary_u,zu,naU,nbU,bU,biU,boundary_d,zd,naD,nbD,bD,biD, \
                eta_u,eta_d, \
                pre_streaming_BC, \
                rank)

                totalTimeCalcs1 += (time.time()-t0)

#               Finalize getting of edges

                t0 = time.time()

                get_edges_end()

                totalTimeEdges2 += (time.time()-t0)
                totalTimeEdges   = totalTimeEdges1 + totalTimeEdges2

#               Calculations in domain edges

                t0 = time.time()

                overlap  = True ; interior = False

                TLBM_fast(na,D,nx_i,ny_i,nz_i,mx_i,my_i,mz_i,xedge,yedge,zedge, \
                w,c,shifted_distribution,overlap,interior, \
                rho,eta,u, \
                f_last,f,g_last,g, \
                x_xa,y_ya,z_za, \
                nu_f, tau_g, tau_h, BETA_G, \
                Temperature, Rheology, Nu_formula, \
                nu_l_min, nu_l_min_cap, nu_l_max_cap, gamma_T, Tau_y, N, Delta_g, \
                nu_f_last, tau_f_last, nu_total, \
                nu_star, \
                P1_thresh, P2_thresh, \
                M1_pd, M2_pd, M3_pd, M4_pd, \
                M1_bs, (np_float(2.0) if bs_capped else M2_bs), M3_bs, M4_bs, \
                rho0_pd, rho_melt_pd, M_ext_pd, \
                rho0_bs, rho_melt_bs, M_ext_bs, \
                M_ext_bs_LP_TTG, M_ext_bs_MP_TTG, M_ext_bs_HP_TTG, \
                mask_bs, mask_pd, \
                g_l, \
                T_nd, \
                cumulative_extracted_pd, cumulative_extracted_bs, \
                cumulative_extracted_LP, cumulative_extracted_MP, cumulative_extracted_HP, \
                P0, rho_M, L_ed,\
                h_last, h, \
                DR_2, T0_p, dT_p, z_nd, \
                rho_solid, rho_eff, \
                boundary_l,xl,naL,nbL,bL,biL,boundary_r,xr,naR,nbR,bR,biR, \
                boundary_b,yb,naB,nbB,bB,biB,boundary_f,yf,naF,nbF,bF,biF, \
                boundary_u,zu,naU,nbU,bU,biU,boundary_d,zd,naD,nbD,bD,biD, \
                eta_u,eta_d, \
                pre_streaming_BC, \
                rank)

                totalTimeCalcs2 += (time.time()-t0)
                totalTimeCalcs   = totalTimeCalcs1 + totalTimeCalcs2

#-------------------------------------------------
#       Python code option
#-------------------------------------------------

        elif code == "python":

#           Streaming step

            t0 = time.time()
            if RKTLBM:
               for a in range(na):
                   f_r_new      = f_r_last[:,:,a,:].reshape(_NX_)[_x_xa_[a]]
                   f_r[:,:,a,:] = f_r_new.reshape(ny_i,nz_i,nx_i)
                   f_b_new      = f_b_last[:,:,a,:].reshape(_NX_)[_x_xa_[a]]
                   f_b[:,:,a,:] = f_b_new.reshape(ny_i,nz_i,nx_i)
                   g_new        = g_last[:,:,a,:].reshape(_NX_)[_x_xa_[a]]
                   g[:,:,a,:]   = g_new.reshape(ny_i,nz_i,nx_i)
                   if rheology == "archean2":
                      h_new         = h_last[:, :, a, :].reshape(_NX_)[_x_xa_[a]]
                      h[:, :, a, :] = h_new.reshape(ny_i, nz_i, nx_i)
            else:
                for a in range(na):
                    f_new      = f_last[:,:,a,:].reshape(_NX_)[_x_xa_[a]]
                    f[:,:,a,:] = f_new.reshape(ny_i,nz_i,nx_i)
                    g_new      = g_last[:,:,a,:].reshape(_NX_)[_x_xa_[a]]
                    g[:,:,a,:] = g_new.reshape(ny_i,nz_i,nx_i)
                    if rheology == "archean2":
                       h_new         = h_last[:, :, a, :].reshape(_NX_)[_x_xa_[a]]
                       h[:, :, a, :] = h_new.reshape(ny_i, nz_i, nx_i)


            totalTimeStream += (time.time()-t0)

#           Start of calculations

            t0 = time.time()

#           Apply BC

            apply_BC()

#           Macroscopic properties

            if shifted_distribution:
                drho = np.sum(f,axis=2)
                rho  = 1 + drho
                if rheology == "archean2": 
                   L_ed = np.sum(h,axis=2)

            else:
                if RKTLBM:
                   if shifted_distribution:
                       drho_r = np.sum(f_r,axis=2)
                       drho_b = np.sum(f_b,axis=2)
                       drho_r = np.where(edge_mask, 0, drho_r)  # Avoid overflows at edges
                       drho_b = np.where(edge_mask, 0, drho_b)
                       rho_r  = rho_r_init + drho_r
                       rho_b  = rho_b_init + drho_b
                   else:
                       rho_r = np.sum(f_r,axis=2)
                       rho_b = np.sum(f_b,axis=2)
                   if rheology == "archean2":
                       L_ed = np.sum(h,axis=2)
                else:
                    rho  = np.sum(f,axis=2)
                    if rheology == "archean2":
                       L_ed = np.sum(h,axis=2)


            if shifted_distribution: drho = np.where(edge_mask,0,drho) # Avoid overflows at edges (non-RKTLBM)
                        
            if RKTLBM:
               rho0  = rho_r/rho1 + rho_b/rho2
               rho   = rho_r + rho_b
            else:
               rho = np.where(edge_mask,1,rho)                       # Avoid overflows at edges

            if RKTLBM: eta = np.sum(g,axis=2)/rho0
            else:      eta = np.sum(g,axis=2)/rho

            eta = np.where(edge_mask,0,eta)                            # Avoid overflows at edges
            
            if RKTLBM: u   = np.einsum('yzax,ad->yzdx',(f_r + f_b),c)
            else:      u   = np.einsum('yzax,ad->yzdx',f,c)

            for d in range(D):
                u[:,:,d,:]  = np.where(edge_mask,0,u[:,:,d,:])          # Avoid overflows at edges
                u[:,:,d,:] /= rho[:,:,:]

#           Do calculations for partial melting

            if RKTLBM and rheology == "archean2":
               # Solidii and Liquidii masks for Peridotite/Basalt 
                T             =(eta/((rho_r/rho1+rho_b/rho2)*DR_2) + T0_p + dT_p/2)
                P             = P0 *( rho_r + rho_b)
                eqn1_sol_bs   = 699.85 - (70400/(P*1000 + 354)) + (77800000/(P*1000 + 354)**2) 
                eqn2_sol_bs   = 661.85 + 3.5*P + 6.2*P**2
                mask1_bs      = P <  P1_thresh
                mask2_bs      =(P >= P1_thresh) & (P < 5)
                T_solidus_bs  = eqn1_sol_bs*mask1_bs + eqn2_sol_bs*mask2_bs
                T_liquidus_bs = 1149.85 + 105*P
                T_solidus_pd  = (-5.104 *P**2 + 133*P + 1120.85)
                T_liquidus_pd = (1799.85 + 114 * P)

                # T/P conditions for TTG (to be used)
                            # Eqns A1 and A2 for LP-MP-TTG
                            # P in Kbar not GPa
                T1_cond     = (T  >  (760 - (60 * (P - 1)**2)) )  & (T < ((1000 - (150 * ((P - 1.2)/ 1.2)**2))))
                P1_cond     = (P  >  (-0.5 *((T - 870)/220))   )  & (P < (1.5   + (0.7 * ((T - 700)/ 200))))
                mask_lpttg  =  P  <  1
                mask_mpttg  =  P  >= 1

                LP_TTG      = T1_cond & P1_cond & mask_lpttg
                MP_TTG      = T1_cond & P1_cond & mask_mpttg

                            # Eqns A3 and A4 for for HP_TTG 
                T2_cond      = (T  > 1000)                                & (T < 1100 + 50 *((P-3.5)/3.5)**2)
                P2_cond      = (P  > (2.35 + 0.15 * ((T-1000)/100)))   & (P < 5)
                HP_TTG       = T2_cond & P2_cond
                if lp_gatev < 0.05: LP_TTG = LP_TTG & False
                if mp_gatev < 0.05: MP_TTG = MP_TTG & False
                if hp_gatev < 0.05: HP_TTG = HP_TTG & False

                # Melt extract calculations
                mask1_bs_melt              = (T > T_solidus_bs) & (T < T_liquidus_bs) & (T_liquidus_bs > T_solidus_bs)
                mask2_bs_melt              =  T > T_liquidus_bs
                M0_bs                      =  mask1_bs_melt * ((T - T_solidus_bs)/(T_liquidus_bs - T_solidus_bs)) + (1 *  mask2_bs_melt)
                mask1_pd_melt              = (T > T_solidus_pd) & (T < T_liquidus_pd) & (T_liquidus_pd > T_solidus_pd)
                mask2_pd_melt              =  T > T_liquidus_pd
                M0_pd                      =  mask1_pd_melt * ((T - T_solidus_pd)/(T_liquidus_pd - T_solidus_pd)) + (1 *  mask2_pd_melt)
                M_pd                       = np.where( M0_pd - cumulative_extracted_pd > 0, M0_pd - cumulative_extracted_pd, 0)
                M_bs                       = np.where( M0_bs - cumulative_extracted_bs > 0, M0_bs - cumulative_extracted_bs, 0)
                L_den                      = np.where( L_ed > 0.1, L_ed, 0)
                M_ext_pd                   = np.where((M_pd > M2_pd) & ((M_pd - M1_pd) < M4_pd) & ((M_pd - M1_pd) > 0), (M_pd - M1_pd) * (1 - L_den) * mask_pd, 0)
                M_ext_bs                   = np.where((M_bs > (np_float(2.0) if bs_capped else M2_bs)) & ((M_bs - M1_bs) < M4_bs) & ((M_bs - M1_bs) > 0), (M_bs - M1_bs) *      L_den  * mask_bs, 0)                
                M_ext_bs_LP_TTG            = M_ext_bs * LP_TTG
                M_ext_bs_MP_TTG            = M_ext_bs * MP_TTG
                M_ext_bs_HP_TTG            = M_ext_bs * HP_TTG
                mask_extracted             = (cumulative_extracted_bs + M_ext_bs < M4_bs)
                cumulative_extracted_bs    = cumulative_extracted_bs  + (M_ext_bs * mask_extracted)
                cumulative_extracted_LP    = cumulative_extracted_LP  + (M_ext_bs_LP_TTG * mask_extracted) 
                cumulative_extracted_MP    = cumulative_extracted_MP  + (M_ext_bs_MP_TTG * mask_extracted)
                cumulative_extracted_HP    = cumulative_extracted_HP  + (M_ext_bs_HP_TTG * mask_extracted)
                cumulative_extracted_bs                 = np.where((cumulative_extracted_bs + M_ext_bs < M4_bs), cumulative_extracted_bs + M_ext_bs, cumulative_extracted_bs)

#           Apply thermal BC. Note that if there is cooling at the Earth's
#           surface due to radiation, or heating/cooling at the CMB,
#           it can be added here (and similarly in the c code)

            if boundary_u: eta[:,zu,:] = eta_u
            if boundary_d: eta[:,zd,:] = eta_d

#-----------------------------------------------------------------------
#           Equilibrium distributions, and forcing term(s) and collision

#           Boussinesq bouyancy term
            if RKTLBM:  
               gi  = beta_g * eta + 3 * g_comp * (rho-rho2)
               g_r = rho_r * gi
               g_b = rho_b * gi

            else:       gi = rho * beta_g * eta

            if RKTLBM and rheology == "archean2": 
                # Density calculation for Peridotite/Basalt
                P_garnet_in                  = 0.014 * T  - 5.4
                P_plag_out                   = 0.024 * T  + 4
                coef_phase_bs                = np.where(P*10 >= P_garnet_in, np.where(P*10 <= P_plag_out, 0.16 * (P*10 - P_garnet_in)/(P_plag_out - P_garnet_in), 0.16), 0)
                coef_res_bs                  = 0.4  * cumulative_extracted_bs
                coef_res_pd                  = np.where( M_pd < M4_pd, -0.04 * cumulative_extracted_pd,0.4*cumulative_extracted_pd)
                coef_phase_pd                = 0
                rho_solid_pd                 = rho0_pd  * (1 + coef_res_pd + coef_phase_pd)
                rho_solid_bs                 = rho0_bs  * (1 + coef_res_bs + coef_phase_bs)
                rho_solid                    =    L_den  * rho_solid_bs  + (1 - L_den) * rho_solid_pd
                rho_eff                      =    L_den  * (rho_solid_bs * (1 - M_ext_bs) + rho_melt_bs * M_ext_bs) + \
                                             (1 - L_den) * (rho_solid_pd * (1 - M_ext_pd) + rho_melt_pd * M_ext_pd)
                gi                          += ((rho_eff-rho_M)/rho_M)*g_l
                rho_b                       += M_ext_bs 
                rho_r                       -= M_ext_bs
                gi                          += ((rho_b-rho_r)/rho_r)


            u2 = 1.5*np.einsum('yzdx,yzdx->yzx',u,u)
            if RKTLBM:
               for a in range(na):
                   cu    = 3*np.einsum('d,yzdx->yzx',np_float(c[a]),u)
                   cu    = w[a] * (1 + cu + 0.5*cu*cu - u2)
                   if shifted_distribution:
                       # Shifted equilibrium: f̃eq_r = feq_r - rho_r_init * C_r[a]
                       # f̃_r_neq = f̃_r - f̃eq_r = (f_r - rho_r_init*C_r) - (feq_r - rho_r_init*C_r) = f_r - feq_r (unchanged)
                       feq_r[:,:,a,:] = rho_r * (C_r[a]+cu) - rho_r_init * C_r[a]
                       feq_b[:,:,a,:] = rho_b * (C_b[a]+cu) - rho_b_init * C_b[a]
                   else:
                       feq_r[:,:,a,:] = rho_r * (C_r[a]+cu)
                       feq_b[:,:,a,:] = rho_b * (C_b[a]+cu)
                   g_eq[:,:,a,:]  = rho0 * eta * (w[a] + cu)
                   f_r_neq        = f_r - feq_r
                   f_b_neq        = f_b - feq_b                
                   g_neq          = g   - (g_eq + Delta_g[a])
                   if rheology == "archean2":
                      h_eq[:, :, a, :]   = L_ed   * w[a] * cu
                      h_neq = h - h_eq
                                     

            else: 
                for a in range(na):
                    cu = 3*np.einsum('d,yzdx->yzx',np_float(c[a]),u)
                    if shifted_distribution:
                        cu = cu + 0.5*cu*cu - u2
                        f_eq[:,:,a,:] = w[a]*(rho * cu + drho)
                        g_eq[:,:,a,:] = eta * (f_eq[:,:,a,:] + w[a])
                    else:
                        cu = 1 + cu + 0.5*cu*cu - u2
                        f_eq[:,:,a,:] = rho * w[a] * cu
                        g_eq[:,:,a,:] = eta * f_eq[:,:,a,:] + Delta_g[a]
                        if rheology == "archean2":
                           h_eq[:, :, a, :]   = L_ed   * w[a] * cu

                    F_b[:,:,a,:] = w[a] * c[a,Z] * gi

                f_neq = f - f_eq
                g_neq = g - g_eq
                if rheology == "archean2":
                   h_neq = h - h_eq


#           Calculate the linear viscosity term nu_linear

            if temperature == "arrhenius":                  # Arrhenius temperature dependence
                if RKTLBM:
                    V_nd = 0    # Assume V_nd = null (Gabriele) since V_nd not specified in Prasanna papers
                    T    = (eta/((rho_r/rho1 + rho_b/rho2)*DR_2) + T0_p + dT_p/2)
                    T_nd = (T - T0_p)/dT_p
                    arg = (act_E + z_nd*V_nd)/(T_nd+1) - act_E/2 
                    nu_linear = (((rho_r/rho1) * nu_f0_r) + ((rho_b/rho2) * nu_f0_b)) * np.exp(arg)
                else:
                    V_nd = 0    # Assume V_nd = null (Gabriele) since V_nd not specified in Prasanna papers
                    T    = (eta/DR_2 + T0_p + dT_p/2)
                    T_nd = (T - T0_p)/dT_p
                    arg = (act_E + z_nd*V_nd)/(T_nd+1) - act_E/2
                    nu_linear = nu_f * np.exp(arg)


            elif gamma_T > 0 and temperature == "reynolds": # Reynolds temperature dependence
                 if RKTLBM:
                    T    = (eta/((rho_r/rho1 + rho_b/rho2)*DR_2) + T0_p + dT_p/2)
                    arg       = -gamma_T*T
                    nu_linear = (((rho_r/rho1) * nu_f0_r) + ((rho_b/rho2) * nu_f0_b)) * np.exp(arg)
                 else:
                    T         = (eta/DR_2 + T0_p + dT_p/2)
                    arg       = -gamma_T*T
                    nu_linear = nu_f * np.exp(arg)

            else:                                           # No temperature dependence
                if RKTLBM: nu_linear = (((rho_r/rho1) * nu_f0_r) + ((rho_b/rho2) * nu_f0_b))
                else:      nu_linear = nu_f

            nu_linear = np.where(nu_linear>nu_l_min,nu_linear,nu_l_min)

#           Calculate the total viscosity allowing for linear and any plastic terms

#           Case of linear rheology

            if rheology == "linear":

#               Set the total viscosity to nu_linear

                nu_total         = nu_linear
                nu_f_last[:,:,:] = nu_total [:,:,:]

#               Calculate the relaxation time

                tau_f             = nu_total * 3 + 0.5

#           Case of plastic rheology: calculate shear strain rate & shear stress

            else:

                S         = np.zeros((D,D,ny_i,nz_i,nx_i),dtype=np_float)
                for i in range(D):
                    for j in range(D):
                        for a in range(1,na):
                            cc = c[a,i] * c[a,j]
                            if cc != 0: S[i,j] += f_neq[:,:,a,:] * cc
                S_factor    = 1.5/(rho*tau_f_last) 
                for i in range(D):
                    for j in range(D):
                        S[i,j] *=  S_factor
                strain_rate = np.zeros((ny_i,nz_i,nx_i),dtype=np_float)
                for i in range(D):
                    for j in range(D):
                        strain_rate += S[i,j] * S[i,j]

                strain_rate = np.sqrt(strain_rate)           # Strain rate

                Tau         = rho * nu_f_last * strain_rate  # Shear stress

#               Calculate nu_plastic for plastic rheologies

                if nu_formula  == "Tau" :
                   nu_plastic = nu_f * np.power(Tau_y/(Tau+1.e-5*Tau_y),N)

                if nu_formula == "strain_rate":
                   nu_plastic = nu_star + np.power(Tau_y/((2 * strain_rate) + (1.e-5*Tau_y)),N)

#               Ensure nu_plastic is at least nu_l_min

                nu_plastic = np.where(nu_plastic>nu_l_min,nu_plastic,nu_l_min)

#               Calculate the total viscosity

                denominator       = 1/nu_linear + 1/nu_plastic
                nu_total          = 1/denominator
                nu_total          = np.where(nu_total>nu_l_min,nu_total,nu_l_min)
                if not RKTLBM: nu_f_last [:,:,:] = nu_total[:,:,:]

#               Calculate the relaxation time

                tau_f             = nu_total * 3 + 0.5
                if not RKTLBM: tau_f_last[:,:,:] = tau_f[:,:,:]

#           Add the collision terms and Boussinesq term using MRT

            if True and not RKTLBM:
                Delta_f = collision_MRT_py()
                for a in range(na):
                    f[:,:,a,:] += (Delta_f[:,:,a,:] + F_b[:,:,a,:])
                    g[:,:,a,:] -= (g_neq[:,:,a,:]/tau_g)
                    if rheology == "archean2":
                       h[:,:,a,:] -= (h_neq[:,:,a,:]/tau_h)

#           Add the collision terms and Boussinesq term using SRT

            if RKTLBM:
               # Collision already applied by C kernel; Python must not repeat it.
               if D == 2:
                  for a in range(na):
                      cm[a] = np.sqrt(c[a,0]**2 + c[a,1]**2)
               elif D == 3:
                  for a in range(na):
                      cm[a] = np.sqrt(c[a,0]**2 + c[a,1]**2 + c[a,2]**2)
               for d in range(D):
                   F[d,xedge:mx_i-xedge]  = 0
                   for a in range(na):
                       for y in range(yedge,my_i-yedge):
                           yy_ya   = y_ya[a,y]
                           for z in range(zedge,mz_i-zedge):
                               zz_za   = z_za[a,z]
                               xs      = slice(xedge, mx_i - xedge)
                               xx_xa   = x_xa[a, xs]
                               F[d,xs] -= c[a,d] * (rho_r[yy_ya,zz_za,xx_xa] - rho_b[yy_ya,zz_za,xx_xa])
               if D == 2:
                  fm[xedge:mx_i-xedge]   = np.sqrt(F[0,xedge:mx_i-xedge]**2 + F[1,xedge:mx_i-xedge]**2)
                  for a in range(na):
                      xs = slice(xedge, mx_i - xedge)
                      ttmp[xs] = c[a,0] * F[0,xs] + c[a,1] * F[1,xs]
                      cos_lambda[a,xs] = ttmp[xs] / (fm[xs] * cm[a] + 1.e-10)
                      cos_phi[a,xs]    = ttmp[xs] / (fm[xs] + 1.e-10)
               elif D == 3:
                    fm[xedge:mx_i-xedge] = np.sqrt(F[0,xedge:mx_i-xedge]**2 + F[1,xedge:mx_i-xedge]**2 + F[2,xedge:mx_i-xedge]**2)
                    for a in range(na):
                        xs = slice(xedge, mx_i - xedge) 
                        ttmp[xs]         = (c[a,0] * F[0,xs] + c[a,1] * F[1,xs] + c[a,2] * F[2,xs])
                        cos_lambda[a,xs] = ttmp[xs] / (fm[xs] * cm[a] + 1.e-10)
                        cos_phi[a,xs]    = ttmp[xs] / (fm[xs] + 1.e-10)
               for d in range(D):
                   FF[d,:,:,:] = 0
                   for y in range(yedge,my_i-yedge):
                       for z in range(zedge,mz_i-zedge):
                           FF[d,y,z,xedge:mx_i-xedge] = F[d,xedge:mx_i-xedge]
               if A > 0:
                  xs = slice(xedge, mx_i - xedge)
                  for a in range(na):
                          f_r[:,:,a,:]    += A * fm[xs] * (w[a] * cos_phi[a,xs]**2 - B[a])
                          f_b[:,:,a,:]    += A * fm[xs] * (w[a] * cos_phi[a,xs]**2 - B[a])
               
               xs = slice(xedge, mx_i - xedge)
               active = fm[xs] > 2e-3  # same threshold you used
               if np.any(active): 
                  denom = np.where(rho == 0, 1.0, rho)
                  ff_sum[:, :, :, xs] = (f_r[:, :, :, xs] + f_b[:, :, :, xs]) / denom[:, :, xs]
                  for a in range(na):
                      feq[:, :, a, xs] = w[a] * rho[:, :, xs]
                  coef = (interface_thickness * rho_b * rho_r) / (denom**2)
                  for a in range(na):
                      lam = np.zeros_like(feq[:, :, a, xs])
                      lam[:, :, active] = cos_lambda[a, xedge:mx_i-xedge][active]
                      tempp = coef[:, :, xs] * lam * feq[:, :, a, xs]
                      f_r[:, :, a, xs] = rho_r[:, :, xs] * ff_sum[:, :, a, xs] + tempp
                      f_b[:, :, a, xs] = rho_b[:, :, xs] * ff_sum[:, :, a, xs] - tempp
            else:
                for a in range(na):
                    f[:,:,a,:] -= (f_neq[:,:,a,:]/tau_f - F_b[:,:,a,:])
                    g[:,:,a,:] -= (g_neq[:,:,a,:]/tau_g)
                    if rheology == "archean2":
                       h[:,:,a,:] -= (h_neq[:,:,a,:]/tau_h)

            if RKTLBM:
               f_r_tmp = f_r ; f_r = f_r_last ; f_r_last = f_r_tmp
               f_b_tmp = f_b ; f_b = f_b_last ; f_b_last = f_b_tmp
               g_tmp   = g   ; g   = g_last   ; g_last   = g_tmp
               h_tmp   = h   ; h   = h_last   ; h_last   = h_tmp
               hd_tmp  = hd  ; hd  = hd_last  ; hd_last  = hd_tmp
               hb_tmp  = hb  ; hb  = hb_last  ; hb_last  = hb_tmp
             
            totalTimeCalcs += (time.time()-t0)

#   If needed, get the maximum mach number

    t0 = time.time()

    new_mach_max = False
#   With variable dt the nominal counter it advances in float increments, and a
#   single rank disagreeing on these comparisons deadlocks the collectives
#   inside this block (observed on M1_hd_smC: 14 ranks in get_max_field, one
#   elsewhere). Decide the entry on rank 0 and broadcast so it is identical by
#   construction on every rank.
    _do_mach = ((it != it_last_mach    and it >= (it_last_mach    + Inc_mach   )) or
                (it != it_last_vary_dt and it >= (it_last_vary_dt + Inc_vary_dt) and it >= start_vary_dt))
    _do_mach = comm.bcast(_do_mach, root=0)
    if _do_mach:

        it_mach_calc  = it

#       Calculate mach max on core 0 and put it onto all cores
        v_rms_l       = np.sqrt(get_mean_field(np.einsum('yzdx,yzdx->yzx',u,u)))
        v_rms_p       = v_rms_l*dx_p/(dt*dt_p)
        u_mag         = np.sqrt(np.einsum('yzdx,yzdx->yzx',u,u))
        u_mag_max     = get_max_field(np.float32(u_mag))
        mach_max      = np.abs(u_mag_max/C_s)
        new_mach_max  = True
        mach_max = mach_max*np.ones(1,dtype=np.float32)
        comm.Bcast(mach_max,root=0) ; mach_max = mach_max[0]

#       Calculate dT_fluct on core 0 and put it onto all cores

        T_max,T_min = get_max_min_field(np.float32(eta/DR_2))
#        print_field_range("Time = %d, T" %(it),eta/DR_2) # oooOOO, This is for checking and can be removed after checks, why is dT_fluct not near unity???
        dT_fluct    = T_max-T_min
        dT_fluct    = dT_fluct*np.ones(1,dtype=np.float32)/dT_max
        comm.Bcast(dT_fluct,root=0) ; dT_fluct = dT_fluct[0]
        if dT_fluct > dT_limit_to_exit or math.isnan(dT_fluct) or math.isinf(dT_fluct): dT_overflow = True

#       Calculate rho_bar on core 0, put it onto all cores, calculate drho_rms and put it onto core 0
        if RKTLBM:
           rho = np.maximum(rho_r/rho1 + rho_b/rho2, 1e-30)   
        rho_bar      = get_mean_field(np.float32(rho))
        rho_bar      = rho_bar*np.ones(1,dtype=np.float32)
        comm.Bcast(rho_bar,root=0) ; rho_bar = rho_bar[0]
        drho2        = (np.float32(rho)-rho_bar)**2
        drho_rms     = np.sqrt(get_mean_field(drho2))
        new_rho_bar  = True
        # Note: RKTLBM conserves rho_r+rho_b=1 locally at every cell, so rho_bar
        # stays at ~1 and normalizing the full f_r/f_b arrays is unnecessary and costly.

#       Set the mach overflow flag

        if mach_max > mach_limit_to_exit or math.isnan(mach_max) or math.isinf(mach_max): mach_overflow = True

#   Print the Mach number to keep track of accuracy and stability
#   and print the density fluctuations, and a warning if these are high

    _do_mach_bk = (it != it_last_mach and it >= (it_last_mach + Inc_mach))
    _do_mach_bk = comm.bcast(_do_mach_bk, root=0)
    if _do_mach_bk:
        it_mach      = max(0,it_last_mach + Inc_mach)
        it_last_mach = it_mach
        if rank == 0:
            if output_units == "physical":
                time_text = timesteps_to_physical_time_text(t,dt_p,output_units)
                print("Time = %d (%d computed) = %s, mach = %4.2g" %(it_mach,t_computed,time_text,mach_max),end="",flush=True)
            else:
                print("Time = %d (%d computed), mach = %4.2g" %(it_mach,t_computed,mach_max),end="",flush=True)
            if not RKTLBM: print(", rho = %4.2f +/- %6.4f" %(rho_bar,drho_rms),end="",flush=True)
            else:          print(", rho_bar = %4.2f +/- %6.4f" %(rho_bar,drho_rms),end="",flush=True)
            print(", T_fluct = %6.4f" %(dT_fluct/dT_max),flush=True)

#   Print diagnostics every 1 My to compare RKTLBM vs non-RKTLBM

    t_My_now  = it * dt_p / sec_my
    t_My_prev = max(0, it - 1) * dt_p / sec_my
    if int(t_My_now) > int(t_My_prev):
        if RKTLBM:
            rho_diag = np.maximum(rho_r/rho1 + rho_b/rho2, 1e-30)
        else:
            rho_diag = rho
        _bd1 = bmk_diagnostics()
        if rank == 0 and _bd1:
            print("t = %s  Vrms_nd=%9.4g  v_rms=%8.3f cm/yr  Nu=%7.3f  T_bar=%7.4f  T_mean=%7.1f C" % (
                timesteps_to_physical_time_text(t,dt_p,output_units), _bd1["Vrms_nd"], _bd1["v_rms_cm"], _bd1["Nu"], _bd1["T_bar"], _bd1["T_mean_C"]), flush=True)

#   Change the print increment to print the Mach number as it gets high
#   and may be starting to blow up and eventually exit if it is blowing up

    _do_vary = ((mach_max > max_mach_limit and new_mach_max and it >= start_vary_dt) or
                (it != it_last_vary_dt and it >= it_last_vary_dt + Inc_vary_dt and it >= start_vary_dt))
    _do_vary = comm.bcast(_do_vary, root=0)
    if _do_vary:

        it_vary_dt      = max(start_vary_dt,it_last_vary_dt + Inc_vary_dt)
        it_last_vary_dt = it_vary_dt

#       If the maximum mach number is above max_mach_limit or
#       below the min_mach_limit, update the physical time step

#       Decrease dt if above the max mach limit

        if rank < nr:
            if mach_max > max_mach_limit:

                factor      = factor_decrease_dt     # Factor to decrease time step
                if dt*factor < factor_min: factor = factor_min/dt   # land ON the floor, never below

                if dt > factor_min*1.0001 and nu_f0 > nu_l_min:

                    if factor*nu_f0 < nu_l_min:
                        factor = nu_l_min/nu_f0 # Limit to never go below nu_l_min

                    dt_last  = dt
                    dt      *= factor           # Set the new time step increment
                    beta_g  *= factor**2        # Update gravity => Boussinesq factor
                    BETA_G  *= factor**2        # gauge profile tracks beta_g
                    if RKTLBM:
                       nu_f_r *= factor
                       nu_f_b *= factor
                       nu_g0  *= factor
                       tau_g0    = nu_g0   * 3 + 0.5             
                       tau_h0    = nu_g0   * Pr/Peclet * 3 + 0.5
                       g_comp *= factor**2         # Compositional gravity scales as beta_g
                       nu_star *= factor            # Lattice viscosity offset in strain_rate nu_formula
                       tau_g  -= 0.5 ; tau_g *= factor ; tau_g += 0.5   # Thermal relaxation array: nu_g part scales with dt
                       tau_h  -= 0.5 ; tau_h *= factor ; tau_h += 0.5
                    else:
                       nu_f    *= factor           # Update viscosity
                       nu_f0   *= factor           # Update reference viscosity (for print)
                       nu_g    *= factor           # Update thermal diffusivity
                       tau_g    = nu_g * 3 + 0.5   # Update relaxation time for g
                    if rheology != "linear":
                        Tau_y *= factor**2      # Update the yield stress

                    if rank == 0:
                        print("------------------------------------------------------------------",flush=True)
                        print("Time = %d, MACH = %5.3g -> exceeds %5.3g = max_mach_limit" %(it_vary_dt,mach_max,max_mach_limit),flush=True)
                        print("Time = %d, MACH = %5.3g -> decrease dt to %6.4f, new nu_l = %8.6g" %(it_vary_dt,mach_max,dt,nu_f0),flush=True)
                        if nu_f0 <= 1.001*nu_l_min:
                            print("Warning: Cannot decrease dt further as nu_f0 = nu_l_min = %8.6g" %(nu_l_min),flush=True)
                        if factor == factor_min:
                            print("Reached factor_max_vary_dt = %5.2 limit to decrease dt" %(factor_max_vary_dt),flush=True)

#           Increase dt if below the min mach limit

            elif mach_max < min_mach_limit:

                factor      = factor_increase_dt     #  Factor to increase time step

                if dt < factor_max:

                    dt_last  = dt
                    dt      *= factor          # Set the new time step increment
                    beta_g  *= factor**2       # Update gravity => Boussinesq factor
                    BETA_G  *= factor**2        # gauge profile tracks beta_g
                    if RKTLBM:
                       nu_f_r *= factor
                       nu_f_b *= factor
                       nu_g0  *= factor
                       tau_g0    = nu_g0   * 3 + 0.5
                       tau_h0    = nu_g0   * Pr/Peclet * 3 + 0.5  
                       g_comp *= factor**2         # Compositional gravity scales as beta_g
                       nu_star *= factor            # Lattice viscosity offset in strain_rate nu_formula
                       tau_g  -= 0.5 ; tau_g *= factor ; tau_g += 0.5   # Thermal relaxation array: nu_g part scales with dt
                       tau_h  -= 0.5 ; tau_h *= factor ; tau_h += 0.5
                    else:
                       nu_f    *= factor           # Update viscosity
                       nu_f0   *= factor           # Update reference viscosity (for print)
                       nu_g    *= factor           # Update thermal diffusivity
                       tau_g    = nu_g * 3 + 0.5   # Update relaxation time for g
                    if rheology != "linear":
                        Tau_y *= factor**2     # Update the yield stress

                    if rank == 0:
                        print("------------------------------------------------------------------",flush=True)
                        print("Time = %d, MACH = %5.3g -> below %5.3g = min_mach_limit" %(it_vary_dt,mach_max,min_mach_limit),flush=True)
                        print("Time = %d, MACH = %5.3g -> increase dt to %6.4f, new nu_l = %8.6g" %(it_vary_dt,mach_max,dt,nu_f0),flush=True)
                        if factor == factor_max_vary_dt:
                            print("Reached factor_max_vary_dt = %5.2 limit to increase dt" %(factor_max_vary_dt),flush=True)

    if mach_overflow:
        if rank == 0:
            printfile = [sys.stderr, sys.stdout]
            for i in range(len(printfile)):
                print("==================================================================",file=printfile[i],flush=True)
                print("ERROR: Time step = %d, Mach # = %4.2f is too high" %(it,mach_max),file=printfile[i],flush=True)
                print("       Exiting from run for run/project = %s/%s" %(RUN,PROJECT),file=printfile[i],flush=True)
                print("==================================================================",file=printfile[i],flush=True)
                print(" ",file=sys.stderr,flush=True)

            print("------------------------------------------------------------------",flush=True)
            print("NOTE: Try changing parameters to get lower Mach numbers such as:",flush=True)
            print()
            print("      (1) Specify a lower physical timestep using time_per_timestep or",flush=True)
            print("          set a lower factor_dt as doing so will decrease nu_l, and",flush=True)
            print("          hence dt_p, by factor_dt, so lower mach numbers",flush=True)
            print("      (2) Try using variable time steps via max_mach_limit",flush=True)
            print("          and Inc_vary_dt which will decrease dt_p every Inc_vary_dt",flush=True)
            print("          timesteps by a factor of factor_decrease_dt if mach exceeds",flush=True)
            print("          max_mach_limit thus allowing mach to be brought down",flush=True)
            print("      (3) Increase resolution (size of model) which leads to, or allows",flush=True)
            print("          via factor_dt, a smaller dt_p so lower mach numbers",flush=True)
            print("      (4) Decrease lattice gravity g_l which is equal to decreasing viscosity",flush=True)
            print("          and hence dt_p by sqrt(g_l_new/g_l_old), so lower mach numbers",flush=True)
            print("      (5) Decrease lattice viscosity which is equal to decreasing",flush=True)
            print("            dt_p by (nu_l_new/nu_l_old) so lower mach numbers",flush=True)
            print("      (6) Increase nu_l_min which sets nu_total >= nu_l_min for T dependent",flush=True)
            print("          and nonlinear rheologies to avoid excessively low nu_total",flush=True)
            print("          which can lead to too low tau_f which can lead to high mach",flush=True)
            print("==================================================================",flush=True)
            print()
        exit()

    if dT_overflow:
        if rank == 0:
            printfile = [sys.stderr, sys.stdout]
            for i in range(len(printfile)):
                print("==================================================================",file=printfile[i],flush=True)
                print("ERROR: Time step = %d, dT overflow, dT_fluct/dT_max = %6.3g > %6.3g, T E [%6.3g,%6.3g]" %(it,dT_fluct/dT_max,dT_limit_to_exit,T_min,T_max),end="",flush=True)
                print("       Exiting from run for run/project = %s/%s" %(RUN,PROJECT),file=printfile[i],flush=True)
                print("==================================================================",file=printfile[i],flush=True)
                print(" ",file=sys.stderr,flush=True)
        exit()

    totalTimeMach += (time.time()-t0)

#--------------------------------------------------------------------------
#   Output
#--------------------------------------------------------------------------

    if output_data and (it == 0 or (it != it_last_data and it >= (it_last_data + Inc_data))):
        it_data       = max(0,it_last_data + Inc_data)
        it_output     = it_data
        it_last_data  = it_data
        _save_fields_now = (int(round(it_data / max(Inc_data,1))) % max(int(fields_every_n_data),1) == 0)
        write_data    = True
    else:
        write_data = False

    if output_plots and (it == 0 or (it != it_last_plots and it >= (it_last_plots + Inc_plots))):
        it_plots      = max(0,it_last_plots + Inc_plots)
        it_output     = it_plots
        it_last_plots = it_plots
        write_plots   = True
    else:
        write_plots = False

    if write_data or write_plots:

        print_filename = False # Option to print the data and plot filenames that are saved

        t0 = time.time()
        if rank == 0:
            now = datetime.now()
            date_string =  now.strftime("%B %d, %Y, %H:%M:%S")
            print("==================================================================",flush=True)
            print("Writing output: %s, run/project = %s/%s" %(date_string,RUN,PROJECT),flush=True)
            print("------------------------------------------------------------------",flush=True)
            if output_units == "lattice":
                print("Time = %d (%d computed), mach = %6.2g" %(it_mach,t_computed,mach_max),flush=True)
            else:
                time_text = timesteps_to_physical_time_text(t,dt_p,output_units)
                print("Time = %d (%d computed) = %s, mach = %6.2g" %(it_mach,t_computed,time_text,mach_max),flush=True)
            print("------------------------------------------------------------------",flush=True)

#       Write and plot temperature

        if fix_scale_T:
            vmin=T0_p ; vmax=T0_p+dT_p
        else:
            vmin = 0  ; vmax = 0
        if plot_T :
           write_and_plot_field(eta/DR_2,"T",r'$\displaystyle T$',"Degrees (C)",r'$^\circ C$', \
               save=True,plot=plot_T,clip=clip_T,gpow=gpow_T,c0=T0_p+dT_p/2,c1=1,vmin=vmin,vmax=vmax)
           #   Write and plot temperature profile
           plot_field_profile(eta/DR_2,"T",r'$\displaystyle T$',c0=T0_p+dT_p/2,c1=1)

#       Write and plot density
#       rho_plot (float64) needed by both the rho plot and Dynamic_viscosity
        if RKTLBM:
           if rheology == "archean2":
              if melt_extraction:
                 rho_plot = (rho_eff + rho_r * rho1_p).astype(np.float64)
              else:
                 rho_plot = (rho_r * rho1_p + rho_b * ((1 - L_ed) * rho2_p + L_ed * rho3_p)).astype(np.float64)
           else:
              rho_plot = (rho_r*rho1_p + rho_b*rho2_p).astype(np.float64)
        else:
           if rheology == "archean2":
              rho_plot = (rho*(rho0_bs * L_ed + rho0_pd * (1-L_ed))).astype(np.float64)
           else:
              rho_plot = (rho * rho_p).astype(np.float64)

        if plot_rho:
           write_and_plot_field(rho_plot,"rho",r'$\displaystyle \rho$',r'kg/m$^3$',r'kg/m$^3$', \
                                save=True,plot=plot_rho,clip=clip_rho,gpow=gpow_rho)


#       Write and plot |u| and/or velocity vectors

        if write_data or (write_plots and (plot_u or plot_vectors)):
            
            U_x0             = get_field(u[:,:,X,:])
            if D == 3: U_y0  = get_field(u[:,:,Y,:])
            U_z0             = get_field(u[:,:,Z,:])

        if write_data and rank == 0:
            u_x_min = np.min(U_x0) ; u_x_max = np.max(U_x0)
            u_z_min = np.min(U_z0) ; u_z_max = np.max(U_z0)
            print("------------------------------------------------------------------",flush=True)
            print("----------------------------- Velocity ---------------------------",flush=True)
            print("Saving u_x E [%9.4g" %(u_x_min),",%9.3g" %(u_x_max),"]",flush=True)
            print("Saving u_z E [%9.4g" %(u_z_min),",%9.3g" %(u_z_max),"]",flush=True)
            if globals().get("_save_fields_now", True):
             filename = data_dir + "u_x_%s" %(timetext(it_data,timescale_Inc_data,dt_p))     ; np.save(filename,U_x0)
            if globals().get("_save_fields_now", True):
             filename = data_dir + "u_z_%s" %(timetext(it_data,timescale_Inc_data,dt_p))     ; np.save(filename,U_z0)
            if D == 3:
                u_y_min = np.min(U_y0) ; u_y_max = np.max(U_y0)
                print("Saving u_y E [%9.4g" %(u_y_min),",%9.3g" %(u_y_max),"]",flush=True)
                if globals().get("_save_fields_now", True):
                 filename = data_dir + "u_y_%s" %(timetext(it_data,timescale_Inc_data,dt_p)) ; np.save(filename,U_y0)

#       Plot |u| and/or velocity vectors

        if rank == 0 and(write_plots and (plot_u or plot_vectors)):

                u_convert, cbar_label_u = get_u_convert(Ra_p,Ra,dx_p,dt_p*dt,output_units=output_units,latex=latex)

                for i in range(len(coordinates)):

                    coordinate = int(coordinates[i])

                    U_x            = np.float32(get_plane(U_x0,plane,coordinate))
                    U_z            = np.float32(get_plane(U_z0,plane,coordinate))
                    if D == 3: U_y = np.float32(get_plane(U_y0,plane,coordinate))

                    U_x            = convert_units(U_x,U_x.dtype,0,u_convert)
                    U_z            = convert_units(U_z,U_z.dtype,0,u_convert)
                    if D == 3: U_y = convert_units(U_y,U_y.dtype,0,u_convert)

                    if   D == 2: FIELD = np.sqrt(U_x**2+U_z**2)
                    elif D == 3: FIELD = np.sqrt(U_x**2+U_y**2+U_z**2)

                    if plot_u:

                        file                = "u"
                        if not latex: title = "|u|"
                        else:         title = r'$\displaystyle |u|$'

                        plot_field(title,file,FIELD,clip_u,gpow_u,zero_mean=False,cbar_label=cbar_label_u,plot_image=True,plot_vectors=vector_overlay)

                    if plot_vectors:

                        file                = "vectors"
                        if not latex: title = "vectors"
                        else:         title = r'$\displaystyle {\bf u}$'

                        plot_field(title,file,FIELD,clip_u,gpow_u,zero_mean=False,cbar_label="",plot_image=False,plot_vectors=True)
                if rank == 0: print("------------------------------------------------------------------",flush=True)

#       RKTLBM output
        if RKTLBM :
           if plot_rho_r or plot_rho_b:
              if rank == 0: print("------------------------------------------------------------------",flush=True)
              if rank == 0: print("------------------------- RKTLBM plots ---------------------------",flush=True)
              write_and_plot_field(rho_r*rho1_p,"rho_r",r'$\displaystyle \rho_{r}$',r'kg/m$^3$',r'kg/m$^3$',save=True,plot=plot_rho_r,plot_vectors=plot_vectors)
              write_and_plot_field(rho_b*rho2_p,"rho_b",r'$\displaystyle \rho_{b}$',r'kg/m$^3$',r'kg/m$^3$',save=True,plot=plot_rho_b,plot_vectors=plot_vectors)
              if rank == 0: print("------------------------------------------------------------------",flush=True)
           if plot_f_r or plot_f_b:
              write_and_plot_field(f_r[:,:,0,:],"f_r",r'$\displaystyle \f_{r} of a = 0$',"",r'',save=True,plot=plot_f_r)
              write_and_plot_field(f_b[:,:,0,:],"f_b",r'$\displaystyle \f_{b} of a = 0$',"",r'',save=True,plot=plot_f_b)
              if rank == 0: print("------------------------------------------------------------------",flush=True)
           if plot_FF_mag or plot_FF_x or plot_FF_z or plot_FF_y:
              if D == 2: F_mag = np.sqrt(FF[X]**2+FF[Z]**2)
              else:      F_mag = np.sqrt(FF[X]**2+FF[Y]**2+FF[Z]**2)
              if plot_FF_mag: write_and_plot_field(F_mag,"F",r'$\displaystyle |{\bf F}|$',"",r'',save=True,plot=plot_FF_mag)
              if plot_FF_x: write_and_plot_field(FF[X],"F_x",r'$\displaystyle F_x$',"",r'',save=True,plot=plot_FF_x)
              if plot_FF_z: write_and_plot_field(FF[Z],"F_z",r'$\displaystyle F_z$',"",r'',save=True,plot=plot_FF_z)

#       Mass diagnostics: per-color and total sums vs t=0 every save interval.
        if RKTLBM:
            BUOY_RHO_R = get_field(rho_r)
            BUOY_RHO_B = get_field(rho_b)
            if rank == 0:
                _sum_r   = float(np.sum(BUOY_RHO_R, dtype=np.float64))
                _sum_b   = float(np.sum(BUOY_RHO_B, dtype=np.float64))
                _sum_tot = _sum_r + _sum_b
                _r0 = getattr(write_and_plot_field, "_buoy_init_r",   _sum_r)
                _b0 = getattr(write_and_plot_field, "_buoy_init_b",   _sum_b)
                _t0 = getattr(write_and_plot_field, "_buoy_init_tot", _sum_tot)
                _dr = _sum_r - _r0
                _db = _sum_b - _b0
                _dt = _sum_tot - _t0
                def _pct_init(cur, init):
                    return f"{cur/init*100:8.4f}% of init" if init > 0 else "       --  (init=0)"
                _rpct_tot = (_sum_r / _sum_tot * 100.0) if _sum_tot > 0 else 0.0
                _bpct_tot = (_sum_b / _sum_tot * 100.0) if _sum_tot > 0 else 0.0
                _tot_drift_pct = (_dt / _t0 * 100.0) if _t0 > 0 else float("nan")
                print("------------------------------------------------------------------",flush=True)
                print(f"-- Mass diagnostics (it={it}):",flush=True)
                print(f"   sum(rho_r) = {_sum_r:+.6e}  delta = {_dr:+.3e}  "
                      f"({_pct_init(_sum_r, _r0)})  ⇒ {_rpct_tot:6.3f}% of total",flush=True)
                print(f"   sum(rho_b) = {_sum_b:+.6e}  delta = {_db:+.3e}  "
                      f"({_pct_init(_sum_b, _b0)})  ⇒ {_bpct_tot:6.3f}% of total",flush=True)
                print(f"   sum(r+b)   = {_sum_tot:+.6e}  delta = {_dt:+.3e}  "
                      f"(total-mass drift = {_tot_drift_pct:+.5f}%)",flush=True)

#       Effective density / solidus density: works for ALL rheologies.
#       For archean2 with melt extraction we use the rho_eff field updated
#       by the C kernel (uses rho_solid + rho_ext_pd + rho_solid_bs split).
#       For everything else we fall back to a simple compositional mix.
        if plot_rho_eff or plot_rho_s:
            if rank == 0: print("------------------------------------------------------------------",flush=True)
            if rank == 0: print("------------------- Effective / solidus density ------------------",flush=True)
            if RKTLBM:
                rho = np.maximum(rho_r + rho_b, 1e-30)
                if rheology == "archean2" and melt_extraction:
                    rho_solid_plot = rho_eff + rho_r * rho1_p
                elif rheology == "archean2":
                    rho_solid_plot = rho_r * rho1_p + rho_b * ((1 - L_ed) * rho2_p + L_ed * rho3_p)
                else:
                    # linear / archean: simple two-color mix using species densities
                    rho_solid_plot = rho_r * rho1_p + rho_b * rho2_p
            else:
                # non-RKTLBM single-fluid: use rho_p as the constant solid density
                rho            = np.maximum(rho_plot, 1e-30) if 'rho_plot' in dir() else rho
                rho_solid_plot = np.full_like(eta, rho_p, dtype=np_float)
            temperature_dependent_density = rho_solid_plot * (1 - beta_p * ((eta/(rho*DR_2) + T0_p+dT_p/2) - T0_p))
            write_and_plot_field(temperature_dependent_density,"rho_eff",r'$\displaystyle \rho_{T}$',r'kg/m$^3$',r'kg/m$^3$',save=True,plot=plot_rho_eff,plot_vectors=plot_vectors)
            write_and_plot_field(rho_solid_plot,"rho_s",r'$\displaystyle \rho_{p}$',r'kg/m$^3$',r'kg/m$^3$',save=True,plot=plot_rho_s,plot_vectors=plot_vectors)

#       Output for archean2 rheology

        if rheology == "archean2":

#           Write and plot the density energy distribution

            if plot_L:
                if rank == 0: print("------------------------------------------------------------------",flush=True)
                write_and_plot_field(L_ed*100,"L",r'$\displaystyle L$',"%",r'Basalt %',save=True,plot=plot_L,plot_vectors=plot_vectors)
#           Write and plot the effective and cummulative melt

            if plot_M_LP or plot_M_MP or plot_M_HP:
                if rank == 0: print("------------------------------------------------------------------",flush=True)
                write_and_plot_field(M_ext_bs_LP_TTG*100,"M_LP",r'$\displaystyle M_{LP}$',"%",r'%',save=True,plot=plot_M_LP)
                write_and_plot_field(M_ext_bs_MP_TTG*100,"M_MP",r'$\displaystyle M_{MP}$',"%",r'%',save=True,plot=plot_M_MP)
                write_and_plot_field(M_ext_bs_HP_TTG*100,"M_HP",r'$\displaystyle M_{HP}$',"%",r'%',save=True,plot=plot_M_HP)
            if plot_ce_pd or plot_ce_bs or plot_de_pd or plot_de_bs:
               if rank == 0: print("------------------------------------------------------------------",flush=True)
               write_and_plot_field(cumulative_extracted_LP*100,"ce_LP",r'$\displaystyle ce_{LP}$',"%",r'%',save=True,plot=plot_ce_LP)
               write_and_plot_field(cumulative_extracted_MP*100,"ce_MP",r'$\displaystyle ce_{MP}$',"%",r'%',save=True,plot=plot_ce_MP)
               write_and_plot_field(cumulative_extracted_HP*100,"ce_HP",r'$\displaystyle ce_{HP}$',"%",r'%',save=True,plot=plot_ce_HP)
               write_and_plot_field(cumulative_extracted_pd*100,"ce_pd",r'$\displaystyle ce_{pd}$',"%",r'%',save=True,plot=plot_ce_pd)
               write_and_plot_field(D_ed*100,"de_pd",r'$\displaystyle D_{ed}$',"%",r'%',save=True,plot=plot_de_pd)
               write_and_plot_field(D_eb*100,"de_bs",r'$\displaystyle D_{eb}$',"%",r'%',save=True,plot=plot_de_bs)
               write_and_plot_field(cumulative_extracted_bs*100,"ce_bs",r'$\displaystyle ce_{bs}$',"%",r'%',save=True,plot=plot_ce_bs)
            if plot_M_bs or plot_M_pd:
               if rank == 0: print("------------------------------------------------------------------",flush=True)
               write_and_plot_field(M_ext_pd*100,"M_ext_pd",r'$\displaystyle M_{pd}$',"%",r'%',save=True,plot=plot_M_pd)
               write_and_plot_field((M_ext_bs_LP_TTG + M_ext_bs_MP_TTG + M_ext_bs_HP_TTG)*100,"M_ext_bs",r'$\displaystyle M_{bs}$',"%",r'%',save=True,plot=plot_M_bs)
  
#           Write and plot T liquidus, T_solidus and Pressure

            if plot_T_liq_pd or plot_T_sol_pd or plot_T_liq_bs or plot_T_sol_bs or plot_P:
                if rank == 0: print("------------------------------------------------------------------",flush=True)
                if plot_T_liq_pd or plot_T_liq_bs:
                    write_and_plot_field(T_liquidus_pd,"T_liq_pd",r'$\displaystyle T_{liq}$',"Degrees (C)",r'$^\circ C$',save=True,plot=plot_T_liq_pd)
                    write_and_plot_field(T_liquidus_bs,"T_liq_bs",r'$\displaystyle T_{liq}$',"Degrees (C)",r'$^\circ C$',save=True,plot=plot_T_liq_bs)
                if plot_T_sol_pd or plot_T_sol_bs:
                    write_and_plot_field(T_solidus_pd,"T_sol_pd",r'$\displaystyle T_{sol}$',"Degrees (C)",r'$^\circ C$',save=True,plot=plot_T_sol_pd)
                    write_and_plot_field(T_solidus_bs,"T_sol_bs",r'$\displaystyle T_{sol}$',"Degrees (C)",r'$^\circ C$',save=True,plot=plot_T_sol_bs)
                if plot_P:
                    write_and_plot_field(P0*(rho_r/rho1 + rho_b/rho2),"P",r'$\displaystyle P$',"Gpa",r'Gpa',save=True,plot=plot_P)

#           Write and plot DI

            if plot_DI:
               Cex_PD = get_field(cumulative_extracted_pd)
               Cex_LP = get_field(cumulative_extracted_LP)
               Cex_MP = get_field(cumulative_extracted_MP)
               Cex_HP = get_field(cumulative_extracted_HP)
               Cex_BS = Cex_LP + Cex_MP + Cex_HP
               RHO_R  = get_field(rho_r)   # must be outside if rank==0 so all ranks participate
               if rank == 0:
                         print("------------------------------------------------------------------",flush=True)
                         print("---------------------- Differentiation Index ---------------------",flush=True)
                         print(f"Time                  = {it * dt_p / sec_my:.4f} My", flush=True)
                         if geometry == "annular":
                            if L_p > 6.6e5:
                               real_area  = (np.abs(theta2-theta1)/360) * (np.pi * 6378**2 - np.pi * 3478**2)
                            else:
                               real_area  = (np.abs(theta2-theta1)/360) * (np.pi * 6378**2 - np.pi * (6378-660)**2)
                            earth_TSA  = 4 * np.pi * (dx_rev[0]/(2*np.pi) * nx)**2
                            total_area = np.sum(dx_rev * dy_rev * nx)
                            total_BS   = np.sum(Cex_BS * dx_rev[np.newaxis, :, np.newaxis] * dy_rev[np.newaxis, :, np.newaxis])
                            total_LP   = np.sum(Cex_LP * dx_rev[np.newaxis, :, np.newaxis] * dy_rev[np.newaxis, :, np.newaxis])
                            total_MP   = np.sum(Cex_MP * dx_rev[np.newaxis, :, np.newaxis] * dy_rev[np.newaxis, :, np.newaxis])
                            total_HP   = np.sum(Cex_HP * dx_rev[np.newaxis, :, np.newaxis] * dy_rev[np.newaxis, :, np.newaxis])
                            total_PD   = np.sum(Cex_PD * dx_rev[np.newaxis, :, np.newaxis] * dy_rev[np.newaxis, :, np.newaxis])
                            DI_bs.append((total_BS/total_area)*100)
                            DI_bs_LP.append((total_LP/total_area)*100)
                            DI_bs_MP.append((total_MP/total_area)*100)
                            DI_bs_HP.append((total_HP/total_area)*100)
                            DI_pd.append((total_PD/total_area)*100)
                            full_circumference = 2 * np.pi * (dx_rev[0]/(2*np.pi) * nx)
                            rho_r_surface      = RHO_R[:, 0, :]                    
                            cont_threshold     = 0.5   # real crust, not a film
                            cont_mask          = (rho_r_surface > cont_threshold).astype(np_float)
                            cell_area_surf     = dx_rev[0] * dy_rev[0]              # km^2 per surface cell
                            area_cont          = np.sum(cont_mask) * cell_area_surf # km^2 of "continent"
                            surf_area_lattice  = rho_r_surface.size * cell_area_surf
                            cont_area_frac     = (area_cont / surf_area_lattice) * 100
                            mass_avg_percent   = np.mean(rho_r_surface, dtype=np.float64) * 100
                            print("Earth's total surface area (code): {:.2e} Km²".format(earth_TSA), flush=True)
                            print("Real area of the Annular Lattice:  {:.2e} Km²".format(real_area), flush=True)
                            print("Total area of the Annular Lattice: {:.2e} Km²".format(total_area), flush=True)
                            print("Area Error:                        {:.2f} %".format((np.abs(total_area-real_area)/real_area)*100, flush=True))
                            print(f"Full circumference (code):         {full_circumference:.2e} Km", flush=True)
                            print(f"Continents total surface area:     {area_cont:.2e} Km²", flush=True)
                            print(f"Continents cover                   {cont_area_frac:.2f} % of zu (cells with rho_r>{cont_threshold:.2f})", flush=True)
                            print(f"Surface mass-average rho_r:        {mass_avg_percent:.2f} %  (TTG mass / surface mass)", flush=True)
                         else: 
                            total_area   = ny * nz * nx
                            if np.sum(Cex_BS) != 0:    DI_bs.append((np.sum(Cex_BS) / total_area)*100)
                            else : DI_bs.append(0)
                            if np.sum(Cex_PD) != 0:    DI_pd.append((np.sum(Cex_PD) / total_area)*100)
                            else : DI_pd.append(0)
                            if np.sum(Cex_LP) != 0: DI_bs_LP.append((np.sum(Cex_LP) / total_area)*100)
                            else:  DI_bs_LP.append(0)
                            if np.sum(Cex_MP) != 0: DI_bs_MP.append((np.sum(Cex_MP) / total_area)*100)
                            else:  DI_bs_MP.append(0)
                            if np.sum(Cex_HP) != 0: DI_bs_HP.append((np.sum(Cex_HP) / total_area)*100)
                            else: DI_bs_HP.append(0)
                         if DI_bs[-1] != 0:
                            Mex_LP_percent.append((DI_bs_LP[-1] / DI_bs[-1]) * 100)
                            Mex_MP_percent.append((DI_bs_MP[-1] / DI_bs[-1]) * 100)
                            Mex_HP_percent.append((DI_bs_HP[-1] / DI_bs[-1]) * 100)
                         else:
                            Mex_LP_percent.append(0)
                            Mex_MP_percent.append(0)
                            Mex_HP_percent.append(0)

                         # Mass-conservation diagnostic: how much of the extracted
                         # basalt mass (cumulative_extracted_bs) actually shows up
                         # as red (rho_r) in the domain. 100% = perfect deposit;
                         # <100% = red lost (deposit cap, BC sink, surface saturation);
                         # >100% = red being manufactured (BC inflation / wetting source).
                         _sum_red  = float(np.sum(RHO_R, dtype=np.float64))   # float16 sum overflows past 65504
                         _sum_cebs = float(np.sum(Cex_BS, dtype=np.float64))
                         _red_pct  = (_sum_red / _sum_cebs * 100.0) if _sum_cebs > 0 else 0.0
                         print(f"Total DI Basalt       = {DI_bs[-1]:.4f} %", flush=True)
                         print(f"red/ce_bs ratio       = {_red_pct:.4f} %  (sum rho_r / sum ce_bs;  ideal=100%)", flush=True)
                         print(f"  sum(rho_r)          = {_sum_red:.6e}", flush=True)
                         print(f"  sum(ce_bs)          = {_sum_cebs:.6e}", flush=True)
                         print(f"LP TTG DI             = {DI_bs_LP[-1]:.4f} %", flush=True)
                         print(f"MP TTG DI             = {DI_bs_MP[-1]:.4f} %", flush=True)
                         print(f"HP TTG DI             = {DI_bs_HP[-1]:.4f} %", flush=True)
                         print(f"LP TTG percent        = {Mex_LP_percent[-1]:.4f} %", flush=True)
                         print(f"MP TTG percent        = {Mex_MP_percent[-1]:.4f} %", flush=True)
                         print(f"HP TTG percent        = {Mex_HP_percent[-1]:.4f} %", flush=True)
                         print(f"Total DI Peridotite   = {DI_pd[-1]:.4f} %", flush=True)

                         # Continental growth time series (Jain/Sizova style verification
                         # data). One line per diagnostic step appended to RUN.growth,
                         # following the .bmk pattern of the W code.
                         growth_file = proj_dir + RUN + ".growth"
                         if not os.path.exists(growth_file):
                             with open(growth_file,"w") as _fg:
                                 _fg.write("# it t_My DI_bs DI_LP DI_MP DI_HP DI_pd cont_area_km2 cont_frac surf_rho_r red_over_cebs (percentages except area)\n")
                         if geometry == "annular":
                             _gr_area = area_cont ; _gr_frac = cont_area_frac ; _gr_mass = mass_avg_percent
                         else:
                             _gr_area = 0.0 ; _gr_frac = 0.0 ; _gr_mass = 0.0
                         with open(growth_file,"a") as _fg:
                             _fg.write("%d %12.6g %9.4f %9.4f %9.4f %9.4f %9.4f %12.5g %9.4f %9.4f %9.4f\n" %(it, it*dt_p/sec_my, DI_bs[-1], DI_bs_LP[-1], DI_bs_MP[-1], DI_bs_HP[-1], DI_pd[-1], _gr_area, _gr_frac, _gr_mass, _red_pct))

                         nt_axis       = np.linspace(0,it*dt_p/sec_my,len(DI_bs))
                         if np_int(it*dt_p/sec_my) == 120 :
                            plt.figure(figsize=(4,2))
                            plt.plot(nt_axis, DI_bs, 'r-', label='Total TTG')
                            plt.plot(nt_axis, DI_bs_LP, 'b-', label='LP TTG')
                            plt.plot(nt_axis, DI_bs_MP, '-', color = "orange",label='MP TTG')
                            plt.plot(nt_axis, DI_bs_HP, 'g-', label='HP TTG')
                            plt.xlabel('Time (Myr)')
                            plt.ylabel('DI (%)')
                            plt.grid(True)
                            plt.tight_layout()
                            plt.legend()
                            plt.savefig(plot_dir + f'DI_{int(it * dt_p / sec_my)}Myr.png', dpi=300)
                            plt.close()
                    
            if rank == 0: print("------------------------------------------------------------------",flush=True)

#       Plot nu_total in physical units

        if plot_nu_kinematic or plot_nu_dynamic:
                if rank == 0: print("----------------- Kinematic and Dynamic Viscosity plots ----------",flush=True)
#               Ensure log(nu_total), the field to be plotted, is not Inf
                FIELD = np.where(np.abs(nu_total)<tiny,tiny,nu_total).astype(np.float64)
                if state == "solid" :
                         Kinematic_viscosity = np.nan_to_num(np.log10((Pr_p/Pr_l) * FIELD * dx_p**2/dt_p), nan=0.0, posinf=40.0, neginf=0.0)
                         rho_plot_safe       = np.where(rho_plot > 0, rho_plot, 1.0)
                         Dynamic_viscosity   = np.nan_to_num(Kinematic_viscosity + np.log10(rho_plot_safe), nan=0.0, posinf=40.0, neginf=0.0)
                         write_and_plot_field(Kinematic_viscosity,"nu_kinematic",r'$\log_{10}(\eta_{Kinematic})$',r'$\eta$ (m$^2$/s)',r'$\eta$ (m$^2$/s)', \
                         save=True,plot=plot_nu_kinematic)
                         comm.Barrier()
                         write_and_plot_field(Dynamic_viscosity,"nu_dynamic",r'$\log_{10}(\eta_{Dynamic})$',r'$\eta$ (Pa$\cdot$s)',r'$\eta$ (Pa$\cdot$s)', \
                         save=True,plot=plot_nu_dynamic)
                else :
                         Kinematic_viscosity = np.nan_to_num(np.log10(np.sqrt(Ra_l/Ra_p) * FIELD * dx_p**2/dt_p), nan=0.0, posinf=40.0, neginf=0.0)
                         rho_plot_safe       = np.where(rho_plot > 0, rho_plot, 1.0)
                         Dynamic_viscosity   = np.nan_to_num(Kinematic_viscosity + np.log10(rho_plot_safe), nan=0.0, posinf=40.0, neginf=0.0)
                         write_and_plot_field(Kinematic_viscosity,"nu_kinematic",r'$\log_{10}(\eta_{Kinematic})$',r'$\eta$ (m$^2$/s)',r'$\eta$ (m$^2$/s)',\
                         save=True,plot=plot_nu_kinematic)
                         comm.Barrier()
                         write_and_plot_field(Dynamic_viscosity,"nu_dynamic",r'$\log_{10}(\eta_{Dynamic})$',r'$\eta$ (Pa$\cdot$s)',r'$\eta$ (Pa$\cdot$s)', \
                         save=True,plot=plot_nu_dynamic)
                if it_plots == 0:
                   NU_K_GLOBAL = get_field(Kinematic_viscosity)   # already log10  (m^2/s)
                   NU_D_GLOBAL = get_field(Dynamic_viscosity)     # already log10  (Pa.s)
                   if rank == 0:
                      nu_kz = np.mean(NU_K_GLOBAL, axis=(0,2))    # avg over y, x
                      nu_dz = np.mean(NU_D_GLOBAL, axis=(0,2))
                      fig, ax1 = plt.subplots(figsize=(4.5,6.0))
                      plt.subplots_adjust(left=0.18,bottom=0.18,top=0.85,right=0.86)

                      # Bottom x-axis: kinematic viscosity (m^2/s)
                      l1, = ax1.plot(nu_kz, z_km, "-",  color="tab:blue",
                                     linewidth=1.2, label=r"$\nu$ kinematic")
                      ax1.invert_yaxis()                          # depth=0 at top
                      ax1.set_xlabel(r"$\log_{10}\;\nu$  (m$^2$/s)",
                                     fontsize=10, color="tab:blue")
                      ax1.set_ylabel("Depth (km)", fontsize=10)
                      ax1.tick_params(axis='x', labelcolor="tab:blue", labelsize=8)
                      ax1.tick_params(axis='y', labelsize=8)
                      ax1.grid(True, which="both", linestyle=":", linewidth=0.5)

                      # Top x-axis: dynamic viscosity (Pa.s)
                      ax2 = ax1.twiny()
                      l2, = ax2.plot(nu_dz, z_km, "--", color="tab:red",
                                     linewidth=1.2, label=r"$\eta$ dynamic")
                      ax2.set_xlabel(r"$\log_{10}\;\eta$  (Pa$\cdot$s)",
                                     fontsize=10, color="tab:red")
                      ax2.tick_params(axis='x', labelcolor="tab:red", labelsize=8)

                      # Place legend BELOW the figure so it never overlaps the
                      # profile lines (especially the long high-viscosity step
                      # at the bottom of the domain).
                      ax1.legend(handles=[l1, l2],
                                 loc="upper center",
                                 bbox_to_anchor=(0.5, -0.10),
                                 ncol=2, fontsize=9, frameon=False)

                      if output_units == "physical":
                         time_title = timesteps_to_physical_time_text(t,dt_p,output_units)
                         ax1.set_title(f"Viscosity profiles  {time_title}",
                                       fontsize=10, pad=22)
                      else:
                         ax1.set_title(f"Viscosity profiles  it={it_plots}",
                                       fontsize=10, pad=22)

                      if t < 0:
                         full_filename = plot_dir + "nu_kinematic_z.png"
                      else:
                         full_filename = plot_dir + "nu_kinematic_z_%s.png" %(timetext(it_plots,timescale_Inc_plots,dt_p))
                      try: os.remove(full_filename)
                      except OSError: pass
                      plt.savefig(full_filename, bbox_inches='tight')
                      plt.close()
                      print(f"Plotting profile  nu_kin E [{nu_kz.min():.3g},{nu_kz.max():.3g}]  "
                            f"nu_dyn E [{nu_dz.min():.3g},{nu_dz.max():.3g}]  to {full_filename}",flush=True)

        if plot_Ra:
                if rank == 0: print("------------------------- Rayleigh number plot -------------------",flush=True)
                FIELD_nu = np.where(np.abs(nu_total)<tiny,tiny,nu_total).astype(np.float64)
                if state == "solid":
                    nu_p_local = (Pr_p/Pr_l) * FIELD_nu * dx_p**2/dt_p   # m^2/s
                else:
                    nu_p_local = np.sqrt(Ra_l/Ra_p) * FIELD_nu * dx_p**2/dt_p
                nu_p_safe   = np.where(nu_p_local > 0, nu_p_local, tiny)
                Ra_local    = (g_p * beta_p * dT_p * L_p**3) / (nu_p_safe * ka_p)
                Ra_log10    = np.nan_to_num(np.log10(Ra_local), nan=0.0, posinf=40.0, neginf=0.0)
                write_and_plot_field(Ra_log10,"Ra",r'$\log_{10}(\mathrm{Ra})$',r'$\log_{10}(\mathrm{Ra})$',r'$\log_{10}(\mathrm{Ra})$',\
                                     save=True,plot=plot_Ra)

        if plot_nu_total:
                if rank == 0: print("------------------------- nu total plots -------------------------",flush=True)
#               Ensure log(nu_total), the field to be plotted, is not Inf
                FIELD = np.where(np.abs(nu_total)<tiny,tiny,nu_total).astype(np.float64)
                FIELD = get_field(np.log10(FIELD))
                if rank == 0:
                   file  = "nu_total"
                   if not latex:
                      title      = "log_10(nu_total)"
                      cbar_label = "log_10(lattice units)"
                   else:
                      title      = r'$\displaystyle \log_{10}(\nu_{total})$'
                      cbar_label = r'$\displaystyle \log_{10}$(lattice units)'
                   FIELD = get_plane(FIELD,plane,coordinate)
                   plot_field(title,file,FIELD,clip_nu_total,gpow_nu_total,zero_mean=False,cbar_label=cbar_label,plot_image=True,plot_vectors=False)
                comm.Barrier()
                if rank == 0: print("------------------------------------------------------------------",flush=True)

#       Plot log_10(Tau) in lattice units

        if plot_Tau and code == "python":
#               Ensure log(Tau), the field to be plotted, is not Inf
                tiny  = np.finfo(Tau.dtype).tiny
                FIELD = np.where(np.abs(Tau)<tiny,tiny,Tau)
                FIELD = get_field(np.log10(FIELD))
                file  = "Tau"
                if not latex:
                    title      = "log_10(Tau)"
                    cbar_label = "log_10(lattice units)"
                else:
                    title      = r'$\displaystyle \log_{10}(Tau)$'
                    cbar_label = r'$\displaystyle \log_{10}$(lattice units)'
                if rank == 0:
                    if globals().get("_save_fields_now", True):
                     filename = data_dir + "Tau_%s" %(timetext(it_data,timescale_Inc_data,dt_p)) ; np.save(filename,FIELD)
                    FIELD = get_plane(FIELD,plane,coordinate)
                    plot_field(title,file,FIELD,clip_Tau,gpow_Tau,zero_mean=False,cbar_label=cbar_label,plot_image=True,plot_vectors=False)

        if rank == 0 and it_output < nt:
            print("==================================================================",flush=True)
            totalTimeOutput += time.time() - t0

#--------------------------------------------------------------------------
#   Print timing if necessary
#--------------------------------------------------------------------------

    if it >= nt or (output_timing and it != it_last_timing and it >= (it_last_timing + Inc_timing)):
        it_timing      = min(nt,max(0,it_last_timing + Inc_timing))
        it_last_timing = it_timing
        if it_timing > 0: print_timing()

#--------------------------------------------------------------------------
#   Exit if at end of run
#--------------------------------------------------------------------------

    if it >= nt:
       if rheology == "archean2":
          if plot_DI : 
                     if rank == 0:
                         print("------------------------------------------------------------------",flush=True)
                         print("---------------------- Differentiation Index ---------------------",flush=True)
                         print(f"Time                  = {it * dt_p / sec_my:.4f} My", flush=True)
                         print(f"Total DI Basalt       = {DI_bs[-1]:.4f} %", flush=True)
                         print(f"LP TTG DI             = {DI_bs_LP[-1]:.4f} %", flush=True)
                         print(f"MP TTG DI             = {DI_bs_MP[-1]:.4f} %", flush=True)
                         print(f"HP TTG DI             = {DI_bs_HP[-1]:.4f} %", flush=True)
                         print(f"Total DI Peridotite   = {DI_pd[-1]:.4f} %", flush=True)
                         print(f"LP TTG percent        = {Mex_LP_percent[-1]:.4f} %", flush=True)
                         print(f"MP TTG percent        = {Mex_MP_percent[-1]:.4f} %", flush=True)
                         print(f"HP TTG percent        = {Mex_HP_percent[-1]:.4f} %", flush=True)
                         print("------------------------------------------------------------------",flush=True)
                         print("-------------------------- Plotting DI ---------------------------",flush=True)
                         nt_axis       = np.linspace(0,it*dt_p/sec_my,len(DI_bs))
                         np.save(data_dir +'DI_bs_Total.npy', DI_bs)
                         np.save(data_dir +'DI_bs_LP.npy', DI_bs_LP)
                         np.save(data_dir +'DI_bs_MP.npy', DI_bs_MP)
                         np.save(data_dir +'DI_bs_HP.npy', DI_bs_HP)
                         np.save(data_dir +'Mex_LP_percent.npy', Mex_LP_percent)
                         np.save(data_dir +'Mex_MP_percent.npy', Mex_MP_percent)
                         np.save(data_dir +'Mex_HP_percent.npy', Mex_HP_percent)
                         np.save(data_dir +'nt.npy', nt_axis)
                         print("-------------------------- DI.npy saved --------------------------",flush=True)
                         print("-------------------------- nt.npy saved --------------------------",flush=True)
                         plt.figure(figsize=(4,2))
                         plt.plot(nt_axis, DI_bs, 'r-', label='Total TTG')
                         plt.plot(nt_axis, DI_bs_LP, 'b-', label='LP TTG')
                         plt.plot(nt_axis, DI_bs_MP, '-', color="orange", label='MP TTG')
                         plt.plot(nt_axis, DI_bs_HP, 'g-', label='HP TTG')
                         plt.xlabel('Time (Myr)')
                         plt.ylabel('DI (%)')
                         plt.grid(True)
                         plt.tight_layout()
                         plt.legend()
                         plt.savefig(plot_dir + f'DI_{int(it * dt_p / sec_my)}Myr.png', dpi=300)
                         plt.close()
                         print("-------------------------- DI plot saved -------------------------",flush=True)
                         plt.figure(figsize=(4,2))
                         plt.plot(nt_axis, Mex_LP_percent, 'b-', label='LP TTG')
                         plt.plot(nt_axis, Mex_MP_percent, '-', color="orange", label='MP TTG')
                         plt.plot(nt_axis, Mex_HP_percent, 'r-', label='HP TTG')
                         plt.xlabel('Time (Myr)')
                         plt.ylabel('Relative proportion (%)')
                         plt.grid(True)
                         plt.tight_layout()
                         plt.legend()
                         plt.savefig(plot_dir + f'RP_{int(it * dt_p / sec_my)}Myr.png', dpi=300)
                         plt.close()
                         print("-------------------------- RP plot saved -------------------------",flush=True)
                         print("------------------------------------------------------------------",flush=True)
                         print("==================================================================",flush=True)

       exit()

#--------------------------------------------------------------------------
#   Update the computed time steps and integral dt = LBM time
#--------------------------------------------------------------------------

    t_computed += 1      # Increment the computed time step
    t          += dt     # Integrate dt to get the LBM time even for variable dt

