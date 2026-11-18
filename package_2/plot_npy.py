#===============================================================================
# Generate plot files (png) from data files (npy)
#===============================================================================

import os
from os.path import exists
import matplotlib
matplotlib.use("Agg") # This is because Agg works without an x server, and avoids problems should the x server not be up
import numpy as np
import matplotlib.pyplot as plt
import time
import sys
from sys import argv
from datetime import datetime
from numpy.polynomial import Polynomial as Poly
import math
from scipy.interpolate import CubicSpline

#===============================================================================
# Function to plot a field profile (ie. mean function of z)
#===============================================================================

def plot_and_write_bmk(field,print_filename=False):

    global Rx, Ry, Rz, nxi, nx_i, nyi, ny_i, nzi, nz_i, D, f, c, write_plots, write_data, write_bmk, size, it, timescale_Inc_data, dt_p, rank
    global proj_dir, PROJECT, RUN

#   Get temperature

    if size > 1:

        T = np.float32(get_field(field))                        # Temperature

#       Get V_rms and maximum magnitude of velocity

        u_l            = np.einsum('yzax,ad->yzdx',f,c)         # Lattice velocity
        u2_l           = np.einsum('yzdx,yzdx->yzx',u_l,u_l)    # Lattice velocity squared
        u_l_max        = get_max_field(np.sqrt(u2_l))           # Maximum magnitude of lattice velocity

        u2_l          *= dx_l**D                                # Lattice velocity magnitude weighted by cell volume
        Vrms_l         = np.sqrt(get_mean_field(u2_l))          # Vrms lattice
        Vrms           = Vrms_l*L_l/ka_l                        # Vrms non-dimensional via multiplication by L/ka per geodynamics benchmark norms

    else:

        T        = np.zeros((ny,nz,nx))
        T[0,:,:] = field

        u_l      = np.zeros((ny,nz,D,nx))                       # Lattice velocity
        if D == 2:
            X = 0 ; Y = 1 ; Z = 1
        else:   
            X = 0 ; Y = 1 ; Z = 2
            filename_ux     = data_dir + "u_x_%s.npy" %(timetext(it,timescale_Inc_data,dt_p))
            filename_uy     = data_dir + "u_y_%s.npy" %(timetext(it,timescale_Inc_data,dt_p))
            filename_uz     = data_dir + "u_z_%s.npy" %(timetext(it,timescale_Inc_data,dt_p))
            for d in range(D):
                if d == X           : u_l[:,:,d,:] = np.load(filename_ux)
                if d == Y and D == 3: u_l[:,:,d,:] = np.load(filename_uy)
                if d == Z           : u_l[:,:,d,:] = np.load(filename_uz)

        u2_l           = np.einsum('yzdx,yzdx->yzx',u_l,u_l)    # Lattice velocity squared
        for z in range(nz): u2_l[:,z,:] = dx[z]**D              # Lattice velocity magnitude weighted by cell volume
        u_l_max        = np.max(np.sqrt(u2_l))                  # Maximum magnitude of lattice velocity
        Vrms_l         = np.sqrt(np.sum(u2_l)/(nx*ny*nz))       # Vrms lattice
        Vrms           = Vrms_l*L_l/ka_l                        # Vrms non-dimensional via multiplication by L/ka per geodynamics benchmark norms
        C_s = 1/np.sqrt(3)                                      # Acoustic wave speed

#   Calculate the benchmark data on core 0 and plot and save it

    if rank == 0:

#       Get the text for physical time and timsteps for files and bmk filename

        time_text_title = timesteps_to_time_text(it,dt_p,output_units)
        time_text_file  = timetext(it,timescale_Inc_data,dt_p)

#       Plot filenames

        filename_T_z         = "T_z_%s.png"    %(time_text_file)
        filename_TBL_T       = "TBL_T_%s.png"  %(time_text_file)
        filename_TBL_Nu      = "TBL_Nu_%s.png" %(time_text_file)
        filename_Nu          = "Nu_%s.png"     %(time_text_file)
        full_filename_T_z    = plot_dir + filename_T_z
        full_filename_TBL_T  = plot_dir + filename_TBL_T
        full_filename_TBL_Nu = plot_dir + filename_TBL_Nu
        full_filename_Nu     = plot_dir + filename_Nu
        time_text_file  = timetext(it,timescale_Inc_data,dt_p)

#       Reynolds Number calculation using u_l_max as per Figure C2 of Mora et al. 2026
#       because: (1) Vrms is saved anyway and can be converted to Re_rms, and
#       (2) Re max scatter reveals variability when there is turbulence, so save Re_max

        Ma            = u_l_max/C_s
        Re            = u_l_max*L_l/nu_l

#       Normalize temperature, get the mean temperature profile Tz, and calculate the volumetric mean, T_bar

        T  = (T-T0_p)/dT_p
        Tz = np.zeros(nz)
        for iz in range(nz): Tz[iz] = np.sum(T[:,iz,:])
        Tz /= (nx*ny)

#       Initialize integral via the trapezoid rule and calculate T_bar = int Tz(z) dz

        offset = 0
        T_bar  = 0.5*( (Tz[0]-offset)*dx[0]**D +    (Tz[nz-1]-offset)*dx[nz-1]**D )
        for iz in range(1,nz-1):         T_bar += ( (Tz[iz  ]-offset)*dx[iz  ]**D )
        T_bar /= L_l
        T_bar += offset

#       Write temperature profile range and save the data

        if size == 1: print(flush=True)

        bmk_file = proj_dir + RUN + ".bmk"

        if write_plots:
            Tz_p = T0_p + dT_p * Tz
            Tz_min     = np.min(Tz_p)
            Tz_max     = np.max(Tz_p)
            if dT_p > 1000:
                Tz_min     = int(Tz_min     + 273.5) - 273
                Tz_max     = int(Tz_max     + 273.5) - 273
                print("        Saving    (npy) T_z E [%5d,%5d]"     %(Tz_min,Tz_max),flush=True)
                print("        Plotting  (png) T_z E [%5d,%5d]"     %(Tz_min,Tz_max),flush=True)
            else:
                print("        Saving    (npy) T_z E [%8.3g,%8.3g]" %(Tz_min,Tz_max),flush=True)
                print("        Plotting  (png) T_z E [%8.3g,%8.3g]" %(Tz_min,Tz_max),flush=True)
        if write_data:
            np.save(filename_Tz,Tz)

#       Calculate normalized distance from the surface and the constants needed
#       to calculate the non-dimensional heat flux in cylindrical coordinates

        r    = radii/L_p
        r1   = r[nz-1]
        r2   = r[0   ]
        ff   = r1/r2
        z    = r2 - r
        iz   = np.linspace(0,1,nz)
        cs_z = CubicSpline(iz,z)

#       Define resolution increase for plots and related variables

        res     = 10
        nz_r    = (nz-1)*res+1
        iz_r    = np.linspace(0,1,nz_r)
        z_r     = cs_z(iz_r)
        r_r     = r2 - z_r

#       Calculate dT/dr by finite differences as a weighted average of forward and backward differences
#       in the interior, and by forward or backward differences at the edges

        dT      = (Tz[1:nz] - Tz[:nz-1])
        d2T     = (dT[1:nz-1] - dT[:nz-2])
        dr      =  r[1:nz] - r[:nz-1]
        dr2     = (dr[1:nz-1] - dr[:nz-2])
        z_fd    = (z[:nz-1]+z[1:])/2
        r_fd    = (r[:nz-1]+r[1:])/2
        dTdr_fd = dT  /dr
        d2Tdr2  = d2T/dr2

        norm_factor = (1+ff)/2

#       Calculate the non-dimensional heat flux and find the maxia adjacent the upper and lower surfaces

        q_fd    = r_fd * np.log(ff) * dTdr_fd
        q_fd   *= norm_factor # oooooooooooo
        q_fd    = q_fd[:nz-1]

        q_u_fd  = q_fd          # q_u_fd = Nu_fd versus depth
        q_d_fd  = np.flip(q_fd) # q_d_fd = Nu_fd versus distance from base

        z_flip  = 1 - np.flip(z )
        Tz_flip =     np.flip(Tz)

#       Get indices at max heat flux by FD

        for iu in range(nz//2):
            if q_fd[iu+1] < q_fd[iu]: break

        for id in range(nz-3,nz//2,-1):
            if q_fd[id] < q_fd[id+1]: break
        id += 1

#       Max heat flux finite difference solition: depth & T and Nu

        Zu      = z   [iu] ; Zd      = z [id] # Use upper case as zu used as index in domain for upper boundary
        Tu      = Tz  [iu] ; Td      = Tz[id]
        Nu_u_fd = q_fd[iu] ; Nu_d_fd = q_fd[id]

#       Set Nu's and TBL widths by finite differences

        Nu_fd   = (Nu_u_fd+Nu_d_fd)/2
        lu_fd   = (1/Nu_u_fd)/2
        ld_fd   = (1/Nu_d_fd)/2
        l_fd    = (lu_fd  +ld_fd  )/2

#       Calculate Nu's and TBL widths by polynomial fit

        frac_l    = 1 # Fraction of TBLW to fit
        order     = 3 # Polynomial order
        n_fit     = 4 # Max number of points to fit

        ku_fit    = int(frac_l*lu_fd*L_l+0.5) ; nu_fit = max(2,ku_fit) ; nu_fit = min(nu_fit,n_fit) ; order_u = min(nu_fit-1,order) ; nu_fit_r = (nu_fit-1)*res + 1
        kd_fit    = int(frac_l*ld_fd*L_l+0.5) ; nd_fit = max(2,kd_fit) ; nd_fit = min(nd_fit,n_fit) ; order_d = min(nd_fit-1,order) ; nd_fit_r = (nd_fit-1)*res + 1

        iu1_fit   = iu     ; iu2_fit   = iu1_fit + nu_fit

        id2_fit   = id + 1 ; id1_fit   = id2_fit - nd_fit

        iu1_fit_r = iu1_fit*res ; iu2_fit_r = iu1_fit_r + nu_fit_r
        id1_fit_r = id1_fit*res ; id2_fit_r = id1_fit_r + nd_fit_r

        ru_fit    = r [iu1_fit:iu2_fit] ; rd_fit = r [id1_fit:id2_fit] ; ru_fit_r  =         r_r[iu1_fit_r:iu2_fit_r]
        Tu_fit    = Tz[iu1_fit:iu2_fit] ; Td_fit = Tz[id1_fit:id2_fit] ; rd_fit_r  = np.flip(r_r[id1_fit_r:id2_fit_r])

        xu_fit    = ru_fit              ; xd_fit = 1 - np.flip(rd_fit)
        yu_fit    = Tu_fit              ; yd_fit =     np.flip(Tu_fit)

#       Scale the points to be fit to within a range x E [0,1] and y E [0,1]

        xu_fit_min = np.min(xu_fit)  ; xu_fit_range = np.max(xu_fit) - xu_fit_min ; xu_fit -= xu_fit_min ; xu_fit /= xu_fit_range
        yu_fit_min = np.min(yu_fit)  ; yu_fit_range = np.max(yu_fit) - yu_fit_min ; yu_fit -= yu_fit_min ; yu_fit /= yu_fit_range
        xd_fit_min = np.min(xd_fit)  ; xd_fit_range = np.max(xd_fit) - xd_fit_min ; xd_fit -= xd_fit_min ; xd_fit /= xd_fit_range
        yd_fit_min = np.min(yd_fit)  ; yd_fit_range = np.max(yd_fit) - yd_fit_min ; yd_fit -= yd_fit_min ; yd_fit /= yd_fit_range

        xu_data   = (z - xu_fit_min)/xu_fit_range ; yu_data = (Tz      - yu_fit_min)/yu_fit_range
        xd_data   = (r - xd_fit_min)/xd_fit_range ; yd_data = (Tz_flip - yd_fit_min)/yd_fit_range

#       Cubic spline of x to obtain x at increased resolution, and from that T for plotting purposes

        ixu_fit   = np.linspace(0,1,nu_fit) ; ixu_fit_r = np.linspace(0,1,nu_fit_r) ; cs_xu_fit = CubicSpline(ixu_fit,xu_fit) ; xu_fit_r  = cs_xu_fit(ixu_fit_r)
        ixd_fit   = np.linspace(0,1,nd_fit) ; ixd_fit_r = np.linspace(0,1,nd_fit_r) ; cs_xd_fit = CubicSpline(ixd_fit,xd_fit) ; xd_fit_r  = cs_xd_fit(ixd_fit_r)

#       Polynomial fit of data and get derivatve

        poly_yu   = Poly.fit(xu_fit,yu_fit,deg=order_u) ; yu_r = poly_yu(xu_fit_r) ; poly_yui = poly_yu.deriv(1) ; yui_fit = poly_yui(xu_fit) ; yui_r = poly_yui(xu_fit_r)
        poly_yd   = Poly.fit(xd_fit,yd_fit,deg=order_d) ; yd_r = poly_yd(xd_fit_r) ; poly_ydi = poly_yd.deriv(1) ; ydi_fit = poly_ydi(xd_fit) ; ydi_r = poly_ydi(xd_fit_r)


#       Convert yi to Nu

#        factor_yi_u   = ru_fit   * np.log(ff) * norm_factor ; yui_fit *= factor_yi_u
#        factor_yi_u_r = ru_fit_r * np.log(ff) * norm_factor ; yui_r   *= factor_yi_u_r
#        factor_yi_d   = rd_fit   * np.log(ff) * norm_factor ; ydi_fit *= factor_yi_d
#        factor_yi_d_r = rd_fit_r * np.log(ff) * norm_factor ; ydi_r   *= factor_yi_d_r

#       Rescale back to normalized temperature and Nu

#        xu_fit *= xu_fit_range ; xu_fit += xu_fit_min ; xu_fit_r *= xu_fit_range ; xu_fit_r += xu_fit_min
#        xd_fit *= xd_fit_range ; xd_fit += xd_fit_min ; xd_fit_r *= xd_fit_range ; xd_fit_r += xd_fit_min
#        yu_fit *= yu_fit_range ; yu_fit -= yu_fit_min ; yu_r     *= yu_fit_range ; yu_r     -= yu_fit_min
#        yd_fit *= yd_fit_range ; yd_fit -= yd_fit_min ; yd_r     *= yd_fit_range ; yd_r     -= yd_fit_min
#        yui_fit *= xu_fit_range/yu_fit_range*norm_factor ; yui_r *= xu_fit_range/yu_fit_range
#        ydi_fit *= xd_fit_range/yd_fit_range*norm_factor ; ydi_r *= xd_fit_range/yd_fit_range

        Nu_u_r = yui_r * norm_factor
        Nu_d_r = ydi_r * norm_factor

        Nu_u_fit = yui_fit[0]
        Nu_d_fit = ydi_fit[0]

#       Set new Nu

        Nu_fit        = (Nu_u_fit+Nu_d_fit)/2

        lu_fit        = 1/Nu_u_fit/2
        ld_fit        = 1/Nu_d_fit/2
        l_fit         = (lu_fit  +ld_fit  )/2

        Nu            = Nu_fit    
        TBLW          = l_fit*L_l

#       Use the old method

        poly_u_old = Poly.fit(ru_fit,Tu_fit,deg=order_u) ; dTudr_fit = poly_u_old.deriv(1)
        poly_d_old = Poly.fit(rd_fit,Td_fit,deg=order_d) ; dTddr_fit = poly_d_old.deriv(1)

        dTudr_fit *= norm_factor # oooooooooooo
        dTddr_fit *= norm_factor # oooooooooooo

        qu_fit_old    = ru_fit  * np.log(ff) * dTudr_fit(ru_fit) ; Nu_u_fit_old = qu_fit_old
        qd_fit_old    = rd_fit  * np.log(ff) * dTddr_fit(rd_fit) ; Nu_d_fit_old = qd_fit_old
        Nu_u_fit_old  = qu_fit_old[0] ; Nu_d_fit_old  = qd_fit_old[nd_fit-1]

#       Set old Nu

        Nu_fit_old    = (Nu_u_fit_old+Nu_d_fit_old)/2

        lu_fit_old    = 1/Nu_u_fit_old/2
        ld_fit_old    = 1/Nu_d_fit_old/2
        l_fit_old     = (lu_fit_old  +ld_fit_old  )/2

        Nu_old        = Nu_fit_old
        TBLW_old      = l_fit_old*L_l

        if write_plots:

#            xmax     =  max(lu_fit,ld_fit*ff)

            ymin   = max(np.min(yu_r),np.min(yd_r))
            ymax   = max(np.max(yu_r),np.max(yd_r))
            yrange = ymax-ymin
            space  = 0.05
            ymin  -= space*yrange
            ymax  += space*yrange

            fig = plt.figure(figsize=(7,4.5))
            fig = plt.subplots_adjust(left=0.12,bottom=0.12,top=0.93,right=0.96)

            xmin = 0
            xmax   = 1.5*max(lu_fit,ld_fit*ff)
            xmax   = min(xmax,0.05,3)

#            plt.xlim(xmin,xmax)
#            plt.xlim(0,3)

            title = r'$T_{TBL}$(%s)' %(time_text_title)
            plt.title(title,fontsize=14)

            plt.ylabel(r'$T$',fontsize=12)
            plt.xlabel(r'Distance from boundary (upper: blue, lower: red)',fontsize=12)

            plt.plot(xu_fit_r,-yu_r*yu_fit_range*0.5,"-" ,color="blue"                                    ,label=r'$T_u^{fit} $',linewidth=3)
            plt.plot(xu_fit  ,-yu_fit*yu_fit_range*0.5,"o" ,color="blue",markersize=6,markerfacecolor='blue',label=r'$T_u^{data}$'            )
#            plt.plot(z/xu_fit_range      ,Tz,"--",color="blue"                                                          ,linewidth=1)
            plt.plot(xd_fit_r,np.flip(yd_r  )*yd_fit_range*0.5+0.5,"-" ,color="red"                                     ,label=r'$T_d^{fit}$' ,linewidth=3)
            plt.plot(xd_fit  ,np.flip(yd_fit)*yd_fit_range*0.5+0.5,"o" ,color="red" ,markersize=6,markerfacecolor='red' ,label=r'$T_d^{data}$'            )
#            plt.plot(z_flip/xd_fit_range  ,Tz_flip,"--",color="red"                                                           ,linewidth=1)

            plt.legend(loc="lower right",fontsize=12)

            plt.savefig(full_filename_TBL_T)
            plt.close()

            fig = plt.figure(figsize=(7,4.5))
            fig = plt.subplots_adjust(left=0.12,bottom=0.12,top=0.93,right=0.96)

            ymin   = min(np.min(yui_r),np.min(ydi_r))
            ymin   = min(ymin,np.min(yui_fit)) ; ymin = min(ymin,np.max(ydi_fit))

            ymax   = max(np.max(yui_r),np.max(ydi_r))
            ymax   = max(ymax,np.max(yui_fit)) ; ymax = max(ymax,np.max(ydi_fit))
            yrange = ymax-ymin
            ymin -= space*yrange
            ymax += space*yrange

            z_fd     = z[:nz-1] - (z[1:nz]-z[0:nz-1])/2
            zu_fd    = z_fd
            zd_fd    = 1 - z_fd - (z[id]-z[nz-1])/2

            zu_r    = z_r   [iu1_fit_r:iu2_fit_r]
            zd_r    = 1-(z_r[id1_fit_r:id2_fit_r])

            xu_fit  +=          Zu
            xd_fit  += (z[nz-1]-Zd)

#            plt.xlim(0   ,xmax)
#            plt.ylim(ymin,ymax)

            plt.xlabel(r'Distance from boundary (upper: blue, lower: red)',fontsize=12)
            title = r'$Nu_{TBL}$(%s)' %(time_text_title)
            plt.title(title,fontsize=14)
            plt.ylabel(r'$Nu$',fontsize=12)

            plt.plot(xu_fit_r,yui_r   ,"-" ,color="blue"                                                 ,label=r'$Nu_u^{fit }$',linewidth=2)
            plt.plot(xu_fit  ,yui_fit ,"o" ,color="blue"      ,markersize=6,markerfacecolor='blue'       ,label=r'$Nu_u^{data}$ '           )
#            plt.plot(zu_fd   ,q_u_fd  ,"--s",color="blue"     ,markersize=6,markerfacecolor='none'       ,label=r'$Nu_u^{FD } $',linewidth=1)

            plt.plot(xd_fit_r,ydi_r  ,"-"   ,color="red"                                          ,label=r'$Nu_d^{fit }$',linewidth=2)
            plt.plot(xd_fit  ,ydi_fit,"o"   ,color="red",markersize=6,markerfacecolor='red',label=r'$Nu_d^{data}$'            )
#            plt.plot(zd_fd   ,q_d_fd  ,"--s",color="red",markersize=6,markerfacecolor='none'      ,label=r'$Nu_u^{FD } $',linewidth=1)

            title = r'$Nu_{TBL}(%s)$' %(time_text_title)
            plt.title(title,fontsize=12)
            plt.legend(loc="upper right",fontsize=12)

            plt.savefig(full_filename_TBL_Nu)
            plt.close()

#           Plot Nu and the temperature profile


#           Tangent lines of temperature at heat flux maxima

            Ttu    = np.array([Tz[iu],0.5],dtype=float)
            Ttd    = np.array([Tz[id],0.5],dtype=float)

#            dTudz  = - dTudr_fit(r[iu])
#            dTddz  = - dTddr_fit(r[id])
#            dTu    = 0.5 - Tz[iu]
#            dTd    = 0.5 - Tz[id]
#            ztu_cut = z[iu] + dTu/dTudz
#            ztd_cut = z[id] + dTd/dTddz

            ztu_cut = lu_fit            +  z[iu]
            ztd_cut = z[id] - ld_fit*ff + (z[id]-z[nz-1])

            ztu  = np.array([z[iu],ztu_cut],dtype=float)
            ztd  = np.array([z[id],ztd_cut],dtype=float)

#            z_fd     = z[:nz-1] + (z[1:nz]-z[0:nz-1])/2
#            zu_fd    = z_fd
#
#            zd_fd    = 1 - z_fd

            zu_r    = z_r   [iu1_fit_r:iu2_fit_r]
            zd_r    = 1-(z_r[id1_fit_r:id2_fit_r])

            ymin = np.min(Nu_u_r)  ; ymax =  np.max(Nu_u_r)

            yrange = ymax - ymin
            space  = yrange*0.05
            
            ymin  -= space
            ymax  += space

            fig = plt.figure(figsize=(7,4.5))
            fig = plt.subplots_adjust(left=0.12,bottom=0.12,top=0.93,right=0.96)

            plt.xlim(xmin,xmax)
            plt.ylim(ymin,ymax)

            plt.plot(zu_r [::res],Nu_u_r[::res],"o"  ,color="blue"      ,markersize=6)
            plt.plot(zu_r        ,Nu_u_r       ,"-"  ,color="blue"                   ,linewidth=2,label=r'$Nu_u^{fit}$')
            plt.plot(zu_fd       ,q_u_fd       ,"--s",color="blue"                   ,linewidth=1,label=r'$Nu_u^{FD }$',markersize=6,markerfacecolor='none')

            plt.plot(zd_r [::res],Nu_d_r [::res],"o" ,color="red",markersize=2)
            plt.plot(zd_r        ,Nu_d_r       ,"-"  ,color="red"             ,linewidth=2,label=r'$Nu_l^{fit}$')
            plt.plot(zd_fd       ,q_d_fd       ,"--s",color="red"             ,linewidth=1,label=r'$Nu_l^{FD }$',markersize=6,markerfacecolor='none')

            plt.xlabel(r'Distance from boundary (upper: blue, lower: red)',fontsize=12)
            plt.ylabel(r'$Nu$',fontsize=12)

            if output_units == "physical": title = r'$Nu$(%s)' %(time_text_title)
            else                         : title = r'$Nu$(%d)' %(it)

            plt.title(title,fontsize=12)
            plt.legend(loc="best",fontsize=12)

            plt.savefig(full_filename_Nu)
            plt.close()

            print("        ======================================",flush=True)
            print("        BENCHMARK SUMMARY"                     ,flush=True)
            print("        ======================================",flush=True)
            print("        iu             = %9d   = index upper q max" %(iu),flush=True)
            print("        id             = %9d   = index lower q max" %(id),flush=True)
            print("        --------------------------------------",flush=True)
            print("        ff = dx_u/dx_d = %9.4g = lower/upper cell size ratio" %(dx[nz-1]/dx[0 ]),flush=True)
            print("        ln(ff)         = %9.4g" %(np.log(ff))  ,flush=True)
            print("        ======================================",flush=True)
            print("        Thermal boundary layer width"          ,flush=True)
            print("        ======================================",flush=True)
            print("        lu_fd          = %9.4g lattice units" %(lu_fd *L_l),flush=True)
            print("        lu_fit_old     = %9.4g lattice units" %(lu_fit_old*L_l),flush=True)
            print("        lu_fit         = %9.4g lattice units" %(lu_fit*L_l),flush=True)
            print("        --------------------------------------"         ,flush=True)
            print("        ld_fd          = %9.4g lattice units" %(ld_fd     *L_l),flush=True)
            print("        ld_fit_old     = %9.4g lattice units" %(ld_fit_old*L_l),flush=True)
            print("        ld_fit         = %9.4g lattice units" %(ld_fit    *L_l),flush=True)
            print("        --------------------------------------"         ,flush=True)
            print("        l_fd           = %9.4g lattice units" %(l_fd      *L_l) ,flush=True)
            print("        l_fit_old      = %9.4g lattice units" %(l_fit_old *L_l) ,flush=True)
            print("        l_fit          = %9.4g lattice units" %(l_fit     *L_l) ,flush=True)
            print("        --------------------------------------"      ,flush=True)
            print("        TBLW           = %9.4g lattice units" %(TBLW),flush=True)
            print("        ======================================",flush=True)
            print("        Nusselt Numbers (upper and lower)"     ,flush=True)
            print("        ======================================",flush=True)
            print("        Nu_u_fd        = %9.4g" %(Nu_u_fd)     ,flush=True)
            print("        Nu_u_fit_old   = %9.4g" %(Nu_u_fit_old),flush=True)
            print("        Nu_u_fit       = %9.4g" %(Nu_u_fit)    ,flush=True)
            print("        --------------------------------------",flush=True)
            print("        Nu_d_fd        = %9.4g" %(Nu_d_fd)     ,flush=True)
            print("        Nu_d_fit_old   = %9.4g" %(Nu_d_fit_old),flush=True)
            print("        Nu_d_fit       = %9.4g" %(Nu_d_fit)    ,flush=True)
            print("        --------------------------------------",flush=True)
            print("        Nu_u/Nu_d(fit) = %9.4g" %(Nu_u_fit/Nu_d_fit),flush=True)
            print("        ======================================",flush=True)
            print("        Benchmark results"                     ,flush=True)
            print("        ======================================",flush=True)
            print("        Ma             = %9.4g" %(Ma)          ,flush=True)
            print("        Re             = %9.4g" %(Re)          ,flush=True)
            print("        --------------------------------------",flush=True)
            print("        T_bar          = %9.4g" %(T_bar)       ,flush=True)
            print("        Vrms           = %9.4g" %(Vrms)        ,flush=True)
            print("        --------------------------------------",flush=True)
            print("        Nu_fd          = %9.4g"   %(Nu_fd     ),flush=True)
            print("        Nu_fit_old     = %9.4g"   %(Nu_fit_old),flush=True)
            print("        --------------------------------------",flush=True)
            print("        Nu_fit         = %9.4g = %6.4f x Nu_fd"   %(Nu_fit ,Nu_fit/Nu_fd     ),flush=True)
            print("        Nu_fit         = %9.4g = %6.4f x Nu_fit_old"   %(Nu_fit ,Nu_fit/Nu_fit_old),flush=True)
            print("        --------------------------------------",flush=True)
            print("        Saved %s"             %(filename_T_z   ),flush=True)
            print("        Saved %s"             %(filename_TBL_T ),flush=True)
            print("        Saved %s"             %(filename_TBL_Nu),flush=True)
            print("        Saved %s (old method)"%(filename_Nu    ),flush=True)
            if print_filename and write_bmk:
                print("        --------------------------------------",flush=True)
                print("        Appended to %s" %(bmk_file),flush=True)
            print("        ======================================",flush=True)

#           Set the title

            Ra_power     = np.log10(Ra)
            Ra_int_power = int(Ra_power)
            Ra_rem_power = Ra_power-Ra_int_power
            Ra_factor    = 10**Ra_rem_power

            Pr_power     = np.log10(Pr)
            Pr_int_power = int(Pr_power)
            Pr_rem_power = Pr_power-Pr_int_power
            Pr_factor    = 10**Pr_rem_power

            Nu_power     = np.log10(np.abs(Nu)+1e-3)
            Nu_int_power = int(Nu_power)
            Nu_rem_power = Nu_power-Nu_int_power
            Nu_factor    = 10**Nu_rem_power

            Vrms_power     = np.log10(Vrms)
            Vrms_int_power = int(Vrms_power)
            Vrms_rem_power = Vrms_power-Vrms_int_power
            Vrms_factor    = 10**Vrms_rem_power

            Re_power     = np.log10(Re)
            Re_int_power = int(Re_power)
            Re_rem_power = Re_power-Re_int_power
            Re_factor    = 10**Re_rem_power

            Ma_power     = np.log10(Ma)
            Ma_int_power = int(Ma_power)
            Ma_rem_power = Ma_power-Ma_int_power
            Ma_factor    = 10**Ma_rem_power

            if output_units == "physical":
                title = r'$T$(%s): ' %(time_text_title)
            else:
                title = r'$T$(%d): '  %(it)

            if abs(Ra_factor - 1) < 0.05: title = title + r'$Ra = 10^{%d}$,' %(Ra_int_power)
            else                        : title = title + r'$Ra = %3.1f \times 10^{%d},$' %(Ra_factor,Ra_int_power)

            if abs(Pr_factor - 1) < 0.05: title = title + r'$Pr = 10^{%d},$' %(Pr_int_power)
            else                        : title = title + r'$Pr = %3.1f \times 10^{%d},$' %(Pr_factor,Pr_int_power)
            title = title + r'$A=%5.2f,$' %(A)

            if   Nu < 10  : title = title + r'$Nu=%4.2f,$' %(Nu)
            elif Nu < 100 : title = title + r'$Nu=%4.1f,$' %(Nu)
            elif Nu < 1000: title = title + r'$Nu=%5.1f,$' %(Nu)
            else          : title = title + r'$Nu=%4.2f \times 10^{%d},$' %(Nu_factor,Nu_int_power)

            if   Vrms < 10  : title = title + r'$Vrms=%4.2f,$' %(Vrms)
            elif Vrms < 100 : title = title + r'$Vrms=%4.1f,$' %(Vrms)
            elif Vrms < 1000: title = title + r'$Vrms=%5.1f,$' %(Vrms)
            else          : title = title + r'$Vrms=%4.2f \times 10^{%d},$' %(Vrms_factor,Vrms_int_power)

            if   Re < 10  : title = title + r'$Re=%4.2f,$' %(Re)
            elif Re < 100 : title = title + r'$Re=%4.1f,$' %(Re)
            elif Re < 1000: title = title + r'$Re=%5.1f,$' %(Re)
            else          : title = title + r'$Re=%4.2f \times 10^{%d},$' %(Re_factor,Re_int_power)

            if   Ma < 10  : title = title + r'$Ma=%4.2f$' %(Ma)
            elif Ma < 100 : title = title + r'$Ma=%4.1f$' %(Ma)
            elif Ma < 1000: title = title + r'$Ma=%5.1f$' %(Ma)
            else          : title = title + r'$Ma=%4.2f \times 10^{%d}$' %(Ma_factor,Ma_int_power)

#           Set the figure size and labels

            fig = plt.figure(figsize=(9,4.5))
            fig = plt.subplots_adjust(left=0.08,bottom=0.1,top=0.93,right=0.97)
            plt.title(title,fontsize=12)

#           Plot temperature, its tangent lines at the heat flux maxima,
#           and the temperature mid value, volumetric mean and legend

            plt.plot   (z[:nz//2+1],Tz[:nz//2+1],"o-" ,color="blue",linewidth=1,markersize=1.5,markerfacecolor='none')
            plt.plot   (ztu,Ttu                 ,"--" ,color="blue",label="Upper tangent",linewidth=1)
            plt.plot   (z[nz//2:  ],Tz[nz//2:]  ,"o-" ,color="red" ,linewidth=1,markersize=1.5,markerfacecolor='none')
            plt.plot   (ztd,Ttd                 ,"--" ,color="red" ,label="Lower tangent",linewidth=1)

            plt.axhline(0.5  ,  linestyle="-" ,color="black",label=r'$T_{mid} = 0.5$'                      ,linewidth=1)
            plt.axhline(T_bar,  linestyle="--",color="black",label=r'$\overline{T} \ \ \ \ = %5.3f$' %(T_bar),linewidth=1)
            plt.legend(loc="upper center",fontsize=12)

#           Plot labels and set plot limits

            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.xlabel("Depth",fontsize=12)
            plt.ylabel("Temperature",fontsize=12)
            plt.tick_params(axis='x',labelsize=10)
            plt.tick_params(axis='y',labelsize=10)

#           Set filename and save plot

            try: os.remove(full_filename_T_z)
            except OSError: pass
            plt.savefig(full_filename_T_z)
            plt.close()

#       Append benchmark data to the benchmark file

        if write_bmk:

            f_bmk = open(bmk_file,"a")
            print("        Appending (bmk) Nu = %8.3g, Vrms = %8.3g Re = %8.3g, Ma = %8.3g" %(Nu,Vrms,Re,Ma),flush=True)
#
#            f_bmk.writelines("%s    %010d %8.4g %8.4g %8.4g %8.4g %8.4f %8.2f\n" %(RUN,nt,dt_p,math.log(Ra_l,10),math.log(Ra_p,10),math.log(Pr_l,10),math.log(Pr_p,10),A))
#            04/08/2026 format
#            f_bmk.writelines("%s %010d %8.4g %8.4g %8.4g %8.4g %8.4f %8.2f %8.4g %8.4f %8.3f %5d %8.4g %8.4g %8.4g %8.4g\n" %(nt_spec,nt,dt_p,math.log(Ra_l,10),math.log(Ra_p,10),math.log(Pr_l,10),math.log(Pr_p,10),A,t,Ma,Nu,T_bar,Re_l,Re_p,Vrms,Vrms_p))
#           21 Aug format

            f_bmk.writelines("%6d %s %8.5g %8.5g %8.5g %8.5g %8.5g %8.5g %8.5g %8.5g %10.5g\n" %(t_bmk,timescale_Inc_bmk,math.log(Ra_l,10),math.log(Pr_l,10),A,Ma,T_bar,Nu,Vrms,Re,TBLW))
            f_bmk.close()

    exit()

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
# Function to plot a field
#===============================================================================
# There is an option to plot an image, and to plot_vectors (velocity vectors), so
# the velocity vectors can overlay an image. If an image is plotted, then the
# color of the vectors is as specified by vector_color in the input file.
# If plot_image = False, then the velocity vectors are black
#---------------------------------------------------------------------------------------

def plot_field(title,what,field,clip,gpow,zero_mean=False,cbar_label="",plot_image=True,plot_vectors=False,vmin=0,vmax=0):

    global vector_color, full_filename, coordinate, description_alignment

    field         = np.float32(field)
    field_min     = np.min(field)
    field_max     = np.max(field)

    if what == "T":
        if dT_p > 1000:

#           Get temperature range to nearest int, even if temperature is in Centigrade which can be up to -273 degrees

            field_min = int(field_min+273.5) - 273
            field_max = int(field_max+273.5) - 273
            print("        Plotting  (png) %s   E [%5d,%5d]"      %(what,field_min,field_max),end="",flush=True)
        else:
            print("        Plotting  (png) %s   E [%8.3g,%8.3g] " %(what,field_min,field_max),end="",flush=True)
    else:
        print("        Plotting  (png) %s E [%8.3g,%8.3g] " %(what,field_min,field_max)      ,end="",flush=True)
    if D == 3: print(", plane = %s, coordinate = %d" %(plane,coordinate),end=""                     ,flush=True)
    print(flush=True)

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
    if geometry == "annulus":
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
        if np.max(thetas)-np.min(thetas) < 0.2*np.pi:
             aspect = 0.6
        if np.max(thetas)-np.min(thetas) < 0.4*np.pi:
             aspect = 0.7
        else:
             aspect = 1
        aspect0  = aspect

    else:

        aspect = nx/ny

#   Calculate the plot size and set space required for the title,
#   axis labels and colorbar if it is to be plotted

    if geometry == "cartesian" or (geometry == "annulus" and plot_image):
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
    elif geometry == "annulus":
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

    if plot_colorbar and plot_image: space_right  = 0.6
    elif plot_vectors and geometry == "annulus":
                                     space_top    = 0.25
                                     space_bottom = 0.1
                                     space_left   = 0.1
                                     space_right  = 0.1
    elif plot_vectors:               space_right  = 0.35
    else:                            space_right  = 0.1

    space_left   = 0.55
    space_bottom = 0.39
    space_top    = 0.24
    if geometry == "annulus":
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
            if plot_size == "vsmall": space_top   -= 0.05
        space_bottom -= 0.45
        if aspect > 1.05:
            space_bottom -= 0.1
        else:
            space_top    -= 0.75
# oooooooooooo Get new graphics from Amen 
    if geometry == "annulus":
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
            space_top    = 0
            if plot_size == "vsmall": space_top   -= 0.05
# oooooooooooo???   space_left   = 0.05
#        space_bottom -= 0.45 # ooo???
#        space_top    -= 0.75 # ooo???

    space_top    = 0
    plot_size0 = 3
    xinch0     = plot_size0

# oooooooooooo Check size is good for all cases

#    plot_size0 = 2 + aspect/2
#    xinch0     = plot_size0

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

# oooooooooooo Check this works in all cases

    space_right += 0.1
    space_top   += 0.08
    yinch       *= 1.03

    xsize = xinch + space_left + space_right
    ysize = yinch + space_top  + space_bottom

    if geometry == "annulus":
        if np.max(thetas-np.min(thetas)) <= 1.01*np.pi:
            aspect_cbar       = 40*ysize/xsize      # Nice aspect ratio for the colorbar
            fraction_colorbar = 0.015               # Nice colorbar size fraction for annulus geometry
        else:
            aspect_cbar       = 50*ysize/xsize      # Nice aspect ratio for the colorbar
            fraction_colorbar = 0.02                # Nice colorbar size fraction for annulus geometry
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
        elif geometry == "annulus":
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
    if geometry != "annulus":
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
            plt.grid(False)
            image = plt.imshow(field,extent=[x_min,x_max,z_min,z_max],aspect=1,cmap=colormap,vmin=vmin,vmax=vmax)

    else:
        ax = fig.add_axes([space_left,space_bottom,1-space_left-space_right,1-space_top-space_bottom],polar=True)
        if plot_image:
            plt.xlim(theta1*np.pi/180,theta2*np.pi/180)
            plt.ylim(0,np.max(radii))
            plt.grid(False)
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

#            print(", overlaying velocity vectors",end="",flush=True)

            plt.quiver(X[::interval,::interval],Y[::interval,::interval],U1[::interval,::interval]*uscale,-U2[::interval,0::interval]*uscale,color=col_vectors,width=0.002*vector_width,scale=length/vector_length,scale_units='inches')
        elif geometry == "annulus":
            theta, radius = np.meshgrid(thetas,radii)
            U1 = (-U_z*np.cos(theta) - U_x*np.sin(theta))
            U2 = (-U_z*np.sin(theta) + U_x*np.cos(theta))

            u_mag = np.max(np.sqrt(U1**2+U2**2))
            uscale   = xinch0*L_p/R_p/(u_mag+2*L_p/R_p*xinch0*tiny)
            length   = (nz-1)/interval

#            print(", overlaying velocity vectors",end="",flush=True)

            ax.quiver(theta[::interval,::interval],radius[::interval,::interval],U1[::interval,::interval]*uscale,U2[::interval,0::interval]*uscale,color=col_vectors,width=0.002*vector_width,scale=30/vector_length)

#    if what == "u": print(flush=True)

#   Set and plot the title

    if run_title == "unspecified": this_title = RUN
    else:                          this_title = r'$%s$' %(run_title)

    this_title = str2str(this_title)

    if output_units == "physical":
        time_text_title = "%s" %(timesteps_to_time_text(it,dt_p,output_units))
        time_text_file  = timetext(it,timescale_Inc_data,dt_p)
    else:
        time_text_title = "%d timesteps"    %(it)
        time_text_file  = timetext(it,"timesteps",dt_p)

    if t < 0:
        title = this_title + ": " + title
    else:
        title = this_title + ": " + title + "(%s)" %(time_text_title)

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
        if   geometry == "annulus":
            cbar = fig.colorbar(image,orientation='vertical',fraction=fraction_colorbar,pad=pad,aspect=aspect_cbar)
        elif geometry == "cartesian":
            cbar = plt.colorbar(image,orientation='vertical',fraction=fraction_colorbar,pad=pad,aspect=aspect_cbar)

        image.figure.axes[1].tick_params(axis="y",labelsize=fontsize_axis_tics) 
        cbar.set_label(cbar_label,fontsize=cbar_label_scale*fontsize_axis_labels)

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
                    plt.xlabel(r'$x$ (km \times 1000)',fontsize=fontsize_axis_labels)
                    plt.ylabel(r'$z$ (km \times 1000)',fontsize=fontsize_axis_labels)
            elif plane == "xy":
                if not latex:
                    plt.xlabel("x (km x 1000)",fontsize=fontsize_axis_labels)
                    plt.ylabel("y (km x 1000)",fontsize=fontsize_axis_labels)
                else:
                    plt.xlabel(r'$x$ (km \times 1000)',fontsize=fontsize_axis_labels)
                    plt.ylabel(r'$y$ (km \times 1000)',fontsize=fontsize_axis_labels)
            elif plane == "yz":
                if not latex:
                    plt.xlabel("y (km x 1000)",fontsize=fontsize_axis_labels)
                    plt.ylabel("z (km x 1000)",fontsize=fontsize_axis_labels)
                else:
                    plt.xlabel(r'$y$ (km \times 1000)',fontsize=fontsize_axis_labels)
                    plt.ylabel(r'$z$ (km \times 1000)',fontsize=fontsize_axis_labels)

#       Axes in km for small models such as the upper mantle

        else:
            if plane == "xz":
                if not latex:
                    plt.xlabel("x (km)",fontsize=fontsize_axis_labels)
                    plt.ylabel("z (km)",fontsize=fontsize_axis_labels)
                else:
                    plt.xlabel(r'$x$ (km)',fontsize=fontsize_axis_labels)
                    plt.ylabel(r'$z$ (km)',fontsize=fontsize_axis_labels)
            elif plane == "xy":
                if not latex:
                    plt.xlabel("x (km x 1000)",fontsize=fontsize_axis_labels)
                    plt.ylabel("y (km x 1000)",fontsize=fontsize_axis_labels)
                else:
                    plt.xlabel(r'$x$ (km)',fontsize=fontsize_axis_labels)
                    plt.ylabel(r'$y$ (km)',fontsize=fontsize_axis_labels)
            elif plane == "yz":
                if not latex:
                    plt.xlabel("y (km x 1000)",fontsize=fontsize_axis_labels)
                    plt.ylabel("z (km x 1000)",fontsize=fontsize_axis_labels)
                else:
                    plt.xlabel(r'$y$ (km)',fontsize=fontsize_axis_labels)
                    plt.ylabel(r'$z$ (km)',fontsize=fontsize_axis_labels)

#   Print the run description on the plot

    if   description_alignment == "left":
        x_description = x_min+0.01*(x_max-x_min)
    elif description_alignment == "center":
        x_description = x_max-0.01*(x_max-x_min)
    elif description_alignment == "right":
        x_description = x_max-0.01*(x_max-x_min)
   
    if plot_image == False: col_description = "k"
    else:                   col_description = description_color
    if description_location == "image":
        if description_box:
            plt.text(x_description,z_min-0.14*(z_min-z_max),description,fontsize=fontsize_title,color=col_description,horizontalalignment=description_alignment,fontweight='bold',bbox=dict(boxstyle='square,pad=0.3',facecolor='white',edgecolor='black',linewidth=0.5))
        else:
            if col_description == "white" or col_description == "w": weight = 'black'
            else:                                                    weight = 'semibold'
            plt.text(x_description,z_min-0.09*(z_min-z_max),description,fontsize=fontsize_title,color=col_description,horizontalalignment=description_alignment,fontweight=weight)
    elif description_location == "frame":
        if aspect <= 1.5: yshift = -0.03
        else:             yshift =  0.0
        if description_box:

            plt.figtext(x_description/(x_max-x_min),0.06+yshift,description,fontsize=fontsize_title,color=col_description,horizontalalignment=description_alignment,fontweight='bold',bbox=dict(boxstyle='square,pad=0.3',facecolor='white',edgecolor='black',linewidth=0.5))
        else:
            if col_description == "white" or col_description == "w": weight = 'black'
            else:                                                    weight = 'semibold'
            plt.figtext(x_description/(x_max-x_min),0.04+yshift,description,fontsize=fontsize_title,color=col_description,horizontalalignment=description_alignment,fontweight=weight)

#   Set the filename and save the figure

    if D == 2:
        what_file = what
    else:
        what_file = what + "_plane=" + plane 
        if   plane == "xz": what_file = what_file + "_y=" + str(coordinate)
        elif plane == "xy": what_file = what_file + "_z=" + str(coordinate)
        elif plane == "yz": what_file = what_file + "_x=" + str(coordinate)

    if t < 0: full_filename = plot_dir + what_file + ".png"
    else    : full_filename = plot_dir + what_file + "_%s.png" %(time_text_file)

    try: os.remove(full_filename)
    except OSError: pass
    plt.savefig(full_filename)
    plt.close()
    if print_filename:
        print(", %s" %(full_filename),flush=True)
    else:
        print()


#===============================================================================
# Function to read a data file and plot a field
#===============================================================================

def read_and_plot_field(what,what_latex,cbar,cbar_latex, \
    clip=1,gpow=1,zero_mean=False,plot_image=True,plot_vectors=False,vmin=0,vmax=0,c0=0,c1=1):

    global t, it, jt, nt, coordinate, plot_bmk, size

    print("===============================================================================",flush=True)
    print("Making \"%s\" plot files ..." %(what),flush=True)
    print("-------------------------------------------------------------------------------",flush=True)
    time.sleep(0.2)
    it , it_float , factor_it = timevalue_to_timesteps(it0,timescale_Inc_data,dt_p,rank)
    for it in range(it0,it1+1,inc):
        jt = it
#        if it > nt: it = nt
#        it , it_float , factor_it = timevalue_to_timesteps(jt,timescale_Inc_data,dt_p,rank)
        t = it

        time_text_print = timesteps_to_time_text(it,dt_p,output_units)
        time_text_file  = timetext(it,timescale_Inc_data,dt_p)

        filename = data_dir + what + "_%s.npy" %(time_text_file)

        if os.path.exists(filename):
            if D == 2:
                print("Time = %d = %s, " %(it,time_text_print),end="",flush=True)
            else:
                print("Time = %d = %s, plane = %s, coordinate = %d , " %(it,time_text_print,plane,coordinate),end="",flush=True)
            FIELD = np.load(filename)
            file  = what
            if not latex:
                title      = what       ; cbar_label = cbar
            else:
                title      = what_latex ; cbar_label = cbar_latex

            for i in range(len(coordinates)):

                coordinate = int(coordinates[i])
                FIELDi     = get_plane(FIELD,plane,coordinate)

                if output_units != "lattice":
                    FIELDi = convert_units(FIELDi,FIELDi.dtype,c0,c1)
                else:
                    if not latex: cbar_label = "Lattice units"
                    else:         cbar_label = r'Lattice units'

                plot_field(title,what,FIELDi,clip=clip,gpow=gpow,zero_mean=zero_mean,cbar_label=cbar_label,plot_image=plot_image,plot_vectors=plot_vectors,vmin=vmin,vmax=vmax)
        else:
            print(    "ERROR:   Missing data file at %s, missing file is: " %(time_text_print),flush=True)
            print(    "         %s" %(filename),flush=True)
            print(    "Exiting making plot files for \"%s\" due to missing data files!" %(file),flush=True)
            break
#        if jt >= nt_int: break
        if it >= nt: break

#-----------------------------------------------------------
# Main
#-----------------------------------------------------------

# Ensure some default variables are set

PLANE       = "unspecified"
COORDINATES = "unspecified"
print_filename = True                # Option to print the data and plot filenames that are saved
Ma_sim               = "default"     # Desired Mach number for simulation   if time_per_timestep = "calculate" else "default"
u_p_geo              = "default"     # Expected u_p for geodynamics problem if time_per_timestep = "calculate" else "default"
                                     # If a value is specified, it represents cm/y for state = solid, and m/s for state = magma
description_alignment = "center"     # Description alignment (left, center or right)

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
elif len(argv) == 4:
    script, PROJECT, RUN, COORDINATES = argv
    VALID = True
elif len(argv) == 5:
    script, PROJECT, RUN, PLANE, COORDINATES = argv
    VALID = True
elif len(argv) == 6:
    script, PROJECT, RUN, it0, it1, inc = argv
    VALID = True
elif len(argv) == 7:
    script, PROJECT, RUN, it0, it1, inc, COORDINATES = argv
    VALID = True
elif len(argv) == 8:
    script, PROJECT, RUN, it0, it1, inc, PLANE, COORDINATES = argv
    VALID = True

if not VALID or HELP == "-h":

#   Self documentation option

    print("============================================================================",file=sys.stderr,flush=True)
    if not VALID:
        print("ERROR: Incorrect usage of plot_npy.py",file=sys.stderr,flush=True)
    else:
        print("plot_npy.py documentation",file=sys.stderr,flush=True)
    print("============================================================================",file=sys.stderr,flush=True)

    print("Usage: python3 plot_npy.py arguments (one to seven arguments as below)",file=sys.stderr,flush=True)
    print(file=sys.stderr,flush=True)
    print("       python3 plot_npy.py -h          (print documentation)",file=sys.stderr,flush=True)
    print("       python3 plot_npy.py PROJECT RUN ",file=sys.stderr,flush=True)
    print("       python3 plot_npy.py PROJECT RUN coordinates",file=sys.stderr,flush=True)
    print("       python3 plot_npy.py PROJECT RUN plane coordinates",file=sys.stderr,flush=True)
    print("       python3 plot_npy.py PROJECT RUN t0 t1 inc",file=sys.stderr,flush=True)
    print("       python3 plot_npy.py PROJECT RUN t0 t1 inc coordinates",file=sys.stderr,flush=True)
    print("       python3 plot_npy.py PROJECT RUN t0 t1 inc plane coordinates",file=sys.stderr,flush=True)
    print("----------------------------------------------------------------------------",file=sys.stderr,flush=True)
    print("Notes: (1) PROJECT and RUN are the project & run names",file=sys.stderr,flush=True)
    print("       (2) plane is the plane to plot (xz, xy or yz)",file=sys.stderr,flush=True)
    print("       (3) coordinates are the coordinates of the normal",file=sys.stderr,flush=True)
    print("           axis to the plane to plot (eg. if plane=xy",file=sys.stderr,flush=True)
    print('           coordinates must be of form "z=coord1,coord2,..."',file=sys.stderr,flush=True)
    print("       (4) t0, t1 and inc are respectively the start time, end time & increment",file=sys.stderr,flush=True)
    print("           in the same units as the data (npy) files generated by the GTXLBM run",file=sys.stderr,flush=True)
    print("       (5) if plane and/or coordinates are unspecified, they are as",file=sys.stderr,flush=True)
    print("           specified in the input file",file=sys.stderr,flush=True)
    print("       (6) If the time range and increment are unspecified, all",file=sys.stderr,flush=True)
    print("           data files will be used to generate plots per the input",file=sys.stderr,flush=True)
    print("           file specifications (ie. t0=0, t1=nt and inc=Inc_data)",file=sys.stderr,flush=True)
    print("----------------------------------------------------------------------------",file=sys.stderr,flush=True)
    print("Examples: python3 plot_npy.py -h",file=sys.stderr,flush=True)
    print("          python3 plot_npy.py examples archean                          (make plots from all data files of the archean run of the examples project)",file=sys.stderr,flush=True)
    print("          python3 plot_npy.py examples archean 100,200,2                (make plots from 100 to 200 My at an increment of 2 My)",file=sys.stderr,flush=True)
    print("          python3 plot_npy.py examples hadean10 0 28800 30              (make plots every 30 minites for 20 days)",file=sys.stderr,flush=True)
    print("          python3 plot_npy.py examples mantle3D xy 10,90                (make plots of the xy plane for z indices 10 & 90)",file=sys.stderr,flush=True)
    print("          python3 plot_npy.py examples mantle3D 2000 3000 10 z=10,50,90 (make plots from 2000 to 3000 My of plane per input specs)",file=sys.stderr,flush=True)
    print("============================================================================",file=sys.stderr,flush=True)
    exit()

#==============================================================

print(flush=True)
print("===============================================================================",flush=True)
print("Making plot files: run = %s, project = %s" %(RUN,PROJECT),flush=True) 
print("===============================================================================",flush=True)
time.sleep(0.5)

# Set directory names

proj_dir  = "./projects/" + PROJECT + "/"
plot_dir  = proj_dir + RUN + "/plots/"
run_dir   = proj_dir + "/" + RUN + "/"
data_dir  = run_dir + "data/"
plot_dir  = run_dir + "plots/"

# Set the strings for the input and output file names

in_file    = proj_dir + RUN + ".in"
out_file   = proj_dir + RUN + ".out"
data_files = data_dir + "*.npy"
plot_files = plot_dir + "*.png"

#==============================================================
# Check the directories exist and if necessary make them

#   Exit if the project input directory and input run file does not exist

if not os.path.exists(in_file):
    if rank == 0:
        print("ERROR: project input file for run does not exist",file=sys.stderr,flush=True)
        print("       You must create it to do a run",file=sys.stderr,flush=True)
    exit()

#-------------------------------------------------------------------------------
# Import necessary functions
#-------------------------------------------------------------------------------

import get_parameters
from get_parameters import get_parameters
from get_parameters import str2str

import utilities
from utilities import sec_ky
from utilities import timevalue_to_timesteps
from utilities import timetext
from utilities import timevalue
from utilities import timevalue_to_secs
from utilities import timesteps_to_time_text
from utilities import get_plane
from utilities import get_lattice_viscosity
from utilities import convert_units
from utilities import get_u_convert
from utilities import get_annular_lattice
from utilities import get_R_and_L

#-------------------------------------------------------------------------------
# Initialize some variables and for no plots
#-------------------------------------------------------------------------------

Ma_sim               = "default"   # Desired Mach number for simulation   if time_per_timestep = "calculate" else "default"
u_p_geo              = "default"   # Expected u_p for geodynamics problem if time_per_timestep = "calculate" else "default"
nu_l_min             = 0.0001      # Minimum lattice viscosity
gamma_T              = 0           # Reynolds dependent viscosity gamma
rank                 = 0           # Rank = 0 as plot_npy runs on one core program
plot_T               = False       # Temperature
plot_u               = False       # Velocity
plot_vectors         = False       # Velocity
plot_rho             = False       # Density
plot_rho_eff         = False       # Plot of rho effective
plot_rho_s           = False       # Plot of rho solidus
plot_M               = False       # Plot of the melt extract
plot_ce_pd           = False       # Plot the cumulative extracted peridotite
plot_T_sol           = False       # Plot of T solidus
plot_T_liq           = False       # Plot of T liquidus
plot_P               = False       # Plot of the pressure
plot_nu_total        = False       # Plot nu_total     (only for code = python)
plot_Tau             = False       # Plot shear stress (only for code = python)
fix_scale_T          = False       # Fix the temperature scale or not
optimize_nu_l        = False       # Optimize nu_l for accuracy: True or False

#-------------------------------------------------------------------------------
# Get the parameters from the input file including plot choices and scales
#-------------------------------------------------------------------------------

get_parameters(in_file)
from get_parameters import *

if   float_bits == 32: np_float = np.float32
elif float_bits == 64: np_float = np.float64

if   edge_bits  == 16: np_edge  = np.float16
elif edge_bits  == 32: np_edge  = np.float32
elif edge_bits  == 64: np_edge  = np.float64

if   data_bits  == 16: np_data = np.float16
elif data_bits  == 32: np_data = np.float32

if ny == 1: D = 2
else      : D = 3

run_title   = str2str(run_title)
description = str2str(description)

if latex: matplotlib.rcParams['text.usetex'] = True
else    : matplotlib.rcParams['text.usetex'] = False

#-------------------------------------------------------------------------------
# Initialize for the desired geometry
#-------------------------------------------------------------------------------

# Convert R and L from names to values if needed and initialize for the specified geometry

R_p_input = R_p
L_p_input = L_p
fouts     = [sys.stderr,sys.stdout]

R_p, L_p = get_R_and_L(R_p_input,L_p_input,fouts,print_R_and_L_info=False,rank=rank) # Get R_p and L_p from descriptive names if needed

if geometry == "cartesian":
    L_l   = nz-1
    dx_p  = L_p/L_l
    dx    = np.ones(nx,dtype=int)
    z_p   = dx_p*np.linspace(0,L_p,nz,dtype=np_float)
    A     = (nx-1)/(nz-1)
    radii = np.linspace(R_p,R_p-L_p,nz,dtype=np_float)
    dx    = np.ones(nz,dtype=np_float)
elif geometry == "annulus":
    radii, thetas, nz, L_l, dx, A = get_annular_lattice(R_p,L_p,nx,theta1,theta2,BC_sides,rank,print_annulus_info=False)
    dx_p = L_p/L_l
    annulus_interval_spec = annulus_interval
    if type(annulus_interval_spec) == str and annulus_interval == "default": annulus_interval = max((nz-1)//24,4)

#-------------------------------------------------------------------------------
# Calculate the lattice viscosity, thermal diffusivity, gravity and
# physical time step, and calculate Ra_p and Pr_p
#-------------------------------------------------------------------------------

# Convert to secs if timescale != unspecified

dt_p_input  = dt_p ; eta_p_input = eta_p; nu_l_input = nu_l ; g_l_input = g_l

if state == "solid": pturb = "none"

if eta_p != "unspecified":
    eta_p_spec = eta_p
    nu_p_spec  = eta_p_spec/rho_p
    Ra_p_spec  = g_p*beta_p*dT_p*L_p**3/(ka_p*nu_p_spec)
    eta_p_calc = rho_p*g_p*beta_p*dT_p*L_p**3/(ka_p*Ra)
else:
    Ra_p_spec  = Ra
    eta_p_spec = eta_p
    eta_p_calc = eta_p_spec

Ra_p, Pr_p, eta_p, nu_p, ka_p, Ra_l, Pr_l, g_l, nu_l, ka_l, dx_p, dt_p, pturb = get_lattice_viscosity(state,pturb,Ra,Pr,g_p,beta_p,dT_p,L_p,eta_p,rho_p,ka_p,g_l,L_l,nu_l,optimize_nu_l,nu_l_min,dt_p_input,Ma_sim,u_p_geo,gamma_T,rheology,geometry,print_lattice_info=False,rank=rank)

#--------------------------------------------------------------
# Convert from input time units to timesteps and get nt, and the plot and data time increments
#--------------------------------------------------------------

error = False
if type(nt       ) == str: nt_int       , timescale_nt        , error = timevalue("nt"       ,nt       ,error)
if type(Inc_plots) == str: Inc_plots_int, timescale_Inc_plots , error = timevalue("Inc_plots",Inc_plots,error)
if type(Inc_data ) == str: Inc_data_int , timescale_Inc_data  , error = timevalue("Inc_data" ,Inc_data ,error)

nt       , nt_float       , factor_nt    = timevalue_to_timesteps(nt_int       ,timescale_nt           ,dt_p,rank)
Inc_plots, Inc_plots_float, factor_plots = timevalue_to_timesteps(Inc_plots_int,timescale_Inc_plots    ,dt_p,rank)
Inc_data , Inc_data_float , factor_data  = timevalue_to_timesteps(Inc_data_int ,timescale_Inc_data     ,dt_p,rank)

#-------------------------------------------------------------------------------
# Set the plane, coordinates and number of dimensions
#-------------------------------------------------------------------------------

if PLANE != "unspecified":
    plane = PLANE
if COORDINATES != "unspecified":
    coordinates = COORDINATES
else:
    coordinates = str(coordinates)

coordinates = coordinates.split(",")

#-------------------------------------------------------------------------------
# Set output_plots = True, and the time range and increment
#-------------------------------------------------------------------------------

output_plots = True
write_plots  = True
write_data   = False
write_bmk    = False

# Set to generate plot files for all data if time range not specified

if len(argv) == 3 or len(argv) == 4 or len(argv) == 5:
    it0 = 0
    it1 = nt+Inc_data
    inc = Inc_data
elif len(argv) == 6 or len(argv) == 7 or len(argv) == 8:
    it0 = int(it0)//Inc_data_int
    it1 = int(it1)//Inc_data_int
    inc = max(int(inc),Inc_data_int)//Inc_data_int

    inc = inc*Inc_data_int
    it0 = it0*Inc_data_int
    it1 = it1*Inc_data_int

nt_run = nt_int/Inc_data_int

print("Run = %s, project = %s" %(RUN,PROJECT),flush=True)
print("-------------------------------------------------------------------------------",flush=True)
print("nt         = %d %s" %(nt_int       ,timescale_nt       ),flush=True)
print("Inc_data   = %d %s" %(Inc_data_int ,timescale_Inc_data ),flush=True)
print("-------------------------------------------------------------------------------",flush=True)
print("Generating plots from npy files per specifications to plot_npy",flush=True)
print("-------------------------------------------------------------------------------",flush=True)
print("Start time = %d %s" %(it0,timescale_Inc_data),flush=True)
print("End   time = %d %s" %(it1,timescale_Inc_data),flush=True)
print("Increment  = %d %s" %(inc,timescale_Inc_data),flush=True)

#-------------------------------------------------------------------------------
# Plot temperature
#-------------------------------------------------------------------------------

if plot_T:
    plot_bmk = True
    size = 1
    if fix_scale_T:
        vmin=T0_p ; vmax=T0_p+dT_p
    else:
        vmin = 0  ; vmax = 0
    read_and_plot_field("T",r'$\displaystyle T$',"Degrees (C)",r'$^\circ C$', \
            clip=clip_T,gpow=gpow_T,vmin=vmin,vmax=vmax,c0=T0_p+dT_p/2,c1=1)

#-------------------------------------------------------------------------------
# Plot velocity
#-------------------------------------------------------------------------------

if plot_u or plot_vectors:
    print("===============================================================================",flush=True)
    print("Making \"u\" and/or \"vector\" plot files ...",flush=True)
    print("-------------------------------------------------------------------------------",flush=True)
    time.sleep(0.2)
    time.sleep(0.2)
    for it in range(it0,it1+1,inc):
        jt = it
# oooaaa
#        if it > nt: it = nt
#        it , it_float , factor_it = timevalue_to_timesteps(jt,timescale_Inc_data,dt_p,rank)
        t  = it

        filename_ux     = data_dir + "u_x_%s.npy" %(timetext(it,timescale_Inc_data,dt_p))
        filename_uz     = data_dir + "u_z_%s.npy" %(timetext(it,timescale_Inc_data,dt_p))
        if D == 3:
            filename_uy = data_dir + "u_y_%s.npy" %(timetext(it,timescale_Inc_data,dt_p))

        if os.path.exists(filename_ux) and os.path.exists(filename_uz): ufiles = True
        else:                                                           ufiles = False
        if D == 3 and not os.path.exists(filename_uy):                  ufiles = False

        time_text_print = timesteps_to_time_text(it,dt_p,output_units)

        if ufiles:
            U_x0            = np.load(filename_ux)
            U_z0            = np.load(filename_uz)
            if D == 3: U_y0 = np.load(filename_uy)

            u_convert, cbar_label_u = get_u_convert(Ra_p,Ra,dx_p,dt_p,output_units=output_units,latex=latex)

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

#               Convert velocity to output units & set the colorbar label

                if plot_u:
                    file                = "u"
                    if not latex: title = "|u|"
                    else:         title = r'$\displaystyle |u|$'

                    time_text = timesteps_to_time_text(it,dt_p,output_units)
                    if D == 2:
                        print("Time = %d = %s, " %(it,time_text),end="",flush=True)
                    else:
                        print("Time = %d = %s, plane = %s, coordinate = %d, " %(it,time_text,plane,coordinate),end="",flush=True)
                    plot_field(title,file,FIELD,clip_u,gpow_u,zero_mean=False,cbar_label=cbar_label_u,plot_image=True,plot_vectors=vector_overlay)

                if plot_vectors:
                    file                = "vectors"
                    if not latex: title = "vectors"
                    else:         title = r'$\displaystyle {\bf u}$'

                    if not plot_u:
                        if D == 2:
                            print("Time = %d = %s, " %(it,time_text),end="",flush=True)
                        else:
                            print("Time = %d = %s, plane = %s, coordinate = %d, " %(it,time_text,plane,coordinate),end="",flush=True)
                    plot_field(title,file,FIELD,clip_u,gpow_u,zero_mean=False,cbar_label="",plot_image=False,plot_vectors=True)
        else:
            print(    "ERROR:   Missing data files(s) at %s, missing files are: " %(time_text_print),flush=True)
            print(    "         %s" %(filename_ux),flush=True)
            if D == 3:
                print("  and/or %s" %(filename_uy),flush=True)
            print(    "  and/or %s" %(filename_uz),flush=True)
            print(    "Exiting making plot files for \"%s\" due to missing data files!" %(file),flush=True)
            break

#        if jt >= nt_int: break
        if it >= nt: break

#-------------------------------------------------------------------------------
# Plot density
#-------------------------------------------------------------------------------

if plot_rho:
    read_and_plot_field("rho",r'$\displaystyle \rho$',"kg/m^3",r'kg/m$^3$', \
            clip=clip_rho,gpow=gpow_rho)

#-------------------------------------------------------------------------------
# Plot effective density
#-------------------------------------------------------------------------------

if plot_rho_eff:
    read_and_plot_field("rho_eff",r'$\displaystyle \rho_{eff}$',"kg/m^3",r'kg/m$^3$')

#-------------------------------------------------------------------------------
# Plot solidus density
#-------------------------------------------------------------------------------

if plot_rho_s:
    read_and_plot_field("rho_s",r'$\displaystyle \rho_{s}$',"kg/m^3",r'kg/m$^3$')

#-------------------------------------------------------------------------------
# Plot effective melt
#-------------------------------------------------------------------------------

if plot_M:
    read_and_plot_field("M_eff",r'$\displaystyle M_{eff}$',"%",r'%')

#-------------------------------------------------------------------------------
# Plot cumulative peridotite melt
#-------------------------------------------------------------------------------

if plot_ce_pd:
    read_and_plot_field("ce_pd",r'$\displaystyle ce_{pd}$',"%",r'%')

#-------------------------------------------------------------------------------
# Plot T liquidus
#-------------------------------------------------------------------------------

if plot_T_liq:
    read_and_plot_field("T_liq",r'$\displaystyle T_{liq}$',"Degrees (C)",r'$^\circ C$')

#-------------------------------------------------------------------------------
# Plot T solidus
#-------------------------------------------------------------------------------

if plot_T_sol:
    read_and_plot_field("T_sol",r'$\displaystyle T_{sol}$',"Degrees (C)",r'$^\circ C$')

#-------------------------------------------------------------------------------
# Plot pressure
#-------------------------------------------------------------------------------

if plot_P:
    read_and_plot_field("P",r'$\displaystyle P$',"Gpa",r'Gpa')

now = datetime.now()
date_string =  now.strftime("%B %d, %Y, %H:%M:%S")
print("===============================================================================",flush=True)
print("Run = %s, project = %s: all plot files generated!" %(RUN,PROJECT),flush=True)
print("Finished making plot files at %s!" %(date_string),flush=True)
print("===============================================================================",flush=True)

