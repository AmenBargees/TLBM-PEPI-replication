
#-------------------------------------------------------------------------------
# Function to calculate, plot & save benchmark results: Nu, Vrms, Re, T_bar, Ma & TBLW
#-------------------------------------------------------------------------------

def plot_and_write_bmk(T_n,T_latex,print_filename=False):

    global write_plots, D, f, c, dx, dx_l, zedge, mzi

#   Get the normalized temperature to core 0

    T_n = get_field(T_n)

#   Get Vrms and maximum magnitude of velocity

    u_l      = np.einsum('yzax,ad->yzdx',f,c)      # Lattice velocity
    u2_l     = np.einsum('yzdx,yzdx->yzx',u_l,u_l) # Lattice velocity squared
    u_l_max  = get_max_field(np.sqrt(u2_l))        # Maximum magnitude lattice velocity

    if rank < nr:
        dvol   = dx_l**2                          # Weights = volume
        if test_option == "D" or test_option == "H" or test_option == "I" or test_option == "J" or test_option == "K" or test_option == "M" or test_option == "O":
            v2     = u2_l                          # Squared lattice velocity
        elif test_option == "P":
            v2     = u2_l/dx_l                     # Squared lattice velocity / dx
        elif test_option == "S":
            v2     = (u2_l/dx_l**2)                # New test
        elif test_option == "T":
            v2     = u2_l/dx_l                     # Squared lattice velocity (u_p propto 1/sqrt(dx_p))
        else:
            v2     = (u2_l*dx_l**2)                # Squared lattice velocity scaled to physical velocity
    else:
        dvol   = np.zeros(1,dtype=np_float)
        v2     = np.zeros(1,dtype=np_float)

    Vrms   = np.sqrt(get_sum_field(dvol*v2)/(get_sum_field(dvol)+1e-30))
    Vrms  *= (L_l/ka_l)
    v2_max = get_max_field(v2)
    v2_max = v2_max*np.ones(1,dtype=np_float)
    comm.Bcast(v2_max,root=0) ; v2_max = v2_max[0]
    v2    /= v2_max # Normalize

    v2    = get_field(v2)    # Get normalized squared velocity

#   Calculate the benchmark data on core 0 and plot and save it

    if rank == 0:

#       Get the text for physical time and timsteps for files

        time_text_title = timesteps_to_time_text(it,dt_p,output_units)
        time_text_file  = timetext(it,timescale_Inc_data,dt_p)

#       Maximum Mach number

#       Calculate Mach Number and Reynolds Number from Ma. Use u_l_max to calculate Ma per Figure C2
#       of Mora et al. 2026 because: (1) Vrms is saved anyway and can be converted to Re_rms,
#       and (2) Re max scatter reveals variability due to turbulence, so save Re_max = u_l_max/C_s

        Ma            = u_l_max/C_s
        Re            = u_l_max*L_l/nu_l

#       Get the mean temperature and velocity profiles Tz & Vz, and calculate the volumetric mean temperature T_bar

        Tz = np.zeros(nz) ; Vz = np.zeros(nz)
        for iz in range(nz):
            Tz[iz] = np.sum(T_n[:,iz,:])
            Vz[iz] = np.sum(v2 [:,iz,:])
        Tz /= (nx*ny)
        Vz  = v2_max*L_l/ka_l*np.sqrt(Vz)/(nx*ny)

#       Initialize integral via the trapezoid rule and calculate T_bar = int Tz(z) dz

        T_bar  = 0.5*(Tz[0]*dx[0]**D + Tz[nz-1]*dx[nz-1]**D)
        for iz in range(1,nz-1): T_bar += (Tz[iz]*dx[iz]**D)

        w_bar  = 0.5*(dx[0]**D + dx[nz-1]**D)
        for iz in range(1,nz-1): w_bar += dx[iz]**D

        T_bar /= w_bar

#       Write temperature profile range and save the data

        if write_plots or write_data:
            Tz_min     = np.min(Tz) ; Tz_max     = np.max(Tz)
            Vz_min     = np.min(Vz) ; Vz_max     = np.max(Vz)
            if dT_p > 1000:
                Tz_min     = int(Tz_min     + 273.5) - 273
                Tz_max     = int(Tz_max     + 273.5) - 273
                print("        Plotting  (png) T_z E [%5d,%5d]"     %(Tz_min,Tz_max),flush=True)
                print("        Saving    (npy) T_z E [%5d,%5d]"     %(Tz_min,Tz_max),flush=True)
            else:
                print("        Plotting  (png) T_z E [%8.3g,%8.3g]" %(Tz_min,Tz_max),flush=True)
                print("        Saving    (npy) T_z E [%8.3g,%8.3g]" %(Tz_min,Tz_max),flush=True)
            print("        Plotting  (png) V_z E [%8.3g,%8.3g]" %(Vz_min,Vz_max),flush=True)
            print("        Saving    (npy) V_z E [%8.3g,%8.3g]" %(Vz_min,Vz_max),flush=True)
        if write_data:
            filename = data_dir + "T_z_%s.npy" %(time_text_file) ; np.save(filename,Tz)
            filename = data_dir + "V_z_%s.npy" %(time_text_file) ; np.save(filename,Vz)

#       Calculate normalized distance from the surface and the constants needed
#       to calculate the non-dimensional heat flux in cylindrical coordinates

        r        = radii/L_p
        L        = 1
        r1       = r[nz-1]
        r2       = r[0   ]
        ff       = r1/r2
        z        = r2-r
        iz       = np.linspace(0,1,nz)
        cs_z     = CubicSpline(iz,z)
        z_flip   = np.flip(z)
        Tz_flip  = np.flip(Tz)
        q_cond   = (L/r2)/np.log(r2/r1)

#       Calculate dT/dr by finite differences as a weighted average of forward and backward differences
#       in the interior, and by forward or backward differences at the edges

        dT_fd   = np.zeros(nz-1)
        dr      = np.zeros(nz-1)
        for iz in range(nz-1):
            dT_fd[iz] = Tz[iz+1] - Tz[iz]
            dr   [iz] = r [iz+1] - r [iz]

        r_fd      = (r[:nz-1]+r[1:])/2
        dTdr_fd   = dT_fd/dr

        z_fd      = (z[:nz-1]+z[1:])/2 ; z_fd_flip = 1 - np.flip(z_fd)

#       Calculate the non-dimensional heat flux and find the maxia adjacent the upper and base surfaces

        q_fd    = - (1/q_cond) * (r_fd/r2) * dTdr_fd ; q_fd_flip = np.flip(q_fd)

        for iu in range(nz//2):
            if q_fd[iu+1] < q_fd[iu]: break

        for id in range(nz-3,nz//2,-1):
            if q_fd[id] < q_fd[id+1]: break

        id += 1

        qu_fd   = q_fd[iu] ; qd_fd  = q_fd[id]
        Tu      = Tz[iu]   ; Td     = Tz[id]

#       Calculate Nu's & TBL widths at q maxima adjacent upper & base surfaces and the means by finite differences

        Nu_u_fd =  qu_fd
        Nu_d_fd =  qd_fd
        Nu_fd   = (Nu_u_fd+Nu_d_fd)/2

        lu_fd   = -(1/dTdr_fd[iu])/2 ; TBLW_u_fd = lu_fd/dx[iu]*L_l
        ld_fd   = -(1/dTdr_fd[id])/2 ; TBLW_d_fd = ld_fd/dx[id]*L_l
        l_fd    = (lu_fd  +  ld_fd )/2 ; TBLW_fd   = (TBLW_u_fd+TBLW_d_fd)/2

#       Set Nu and TBLW to fd value for printing to out file

        Nu      = Nu_fd
        TBLW    = TBLW_fd

        if write_plots:

#-------------------------------------------------------------------------------
#           Plot temperature versus depth, its tangent lines at the heat flux
#           maxima, and the temperature mid value and volumetric mean
#-------------------------------------------------------------------------------

#           Set the figure size and labels

            fig = plt.figure(figsize=(9,4.5)) ; fig = plt.subplots_adjust(left=0.07,bottom=0.11,top=0.93,right=0.96)

#           Set upper and base temperature tangent line start and end points

            ztu_cut = z[iu]-(0.5-Tz[iu])/dTdr_fd[iu] ; ztd_cut = z[id] - (0.5-Tz[id])/dTdr_fd[id]

            ztu     = np.array([z [iu],ztu_cut],dtype=float) ; ztd = np.array([z [id],ztd_cut],dtype=float)
            Ttu     = np.array([Tz[iu],0.5    ],dtype=float) ; Ttd = np.array([Tz[id],0.5    ],dtype=float)

#           Title with time, Ra, Pr, A and benchmark results

            if output_units == "physical": title = r'$T$(%s): ' %(time_text_title)
            else                         : title = r'$T$(%d): '  %(it)

            Ra_power = np.log10(Ra) ; Ra_int_power = int(Ra_power) ; Ra_rem_power = Ra_power-Ra_int_power ; Ra_factor    = 10**Ra_rem_power
            Pr_power = np.log10(Pr) ; Pr_int_power = int(Pr_power) ; Pr_rem_power = Pr_power-Pr_int_power ; Pr_factor    = 10**Pr_rem_power

            if abs(Ra_factor - 1) < 0.05: title = title + r'$Ra = 10^{%d}$,' %(Ra_int_power)
            else                        : title = title + r'$Ra = %3.1f \times 10^{%d},$' %(Ra_factor,Ra_int_power)

            if abs(Pr_factor - 1) < 0.05: title = title + r'$Pr = 10^{%d},$' %(Pr_int_power)
            else                        : title = title + r'$Pr = %3.1f \times 10^{%d},$' %(Pr_factor,Pr_int_power)
            title = title + r'$A=%5.2f,$' %(A)

            title = title + get_text_from_float(r'$Nu=$',Nu) + "," + get_text_from_float(r'$V_{rms}=$',Vrms) + "," + get_text_from_float(r'$Re=$',Re) + "," + get_text_from_float(r'$Ma=$',Ma)

            plt.title(title,fontsize=12)

            plt.plot   (z[:nz//2+1],Tz[:nz//2+1],"-" ,color="blue"                                                  ,linewidth=1.5)
            plt.plot   (ztu,Ttu                 ,"--",color="blue",label="Upper tangent"                            ,linewidth=1.5)
            plt.plot   (z[nz//2:  ],Tz[nz//2:]  ,"-" ,color="red"                                                   ,linewidth=1.5)
            plt.plot   (ztd,Ttd                 ,"--",color="red" ,label="Lower tangent"                            ,linewidth=1.5)
            plt.axhline(0.5           ,linestyle="-" ,color="black",label=r'$T_{mid}=0.5$'                          ,linewidth=0.75)
            plt.axhline(T_bar         ,linestyle="--",color="black",label=r'$\overline{T} \ \ \ \ = %5.3f$' %(T_bar),linewidth=0.75)

            plt.legend(loc="upper center",fontsize=12)

#           Axis labels and plot limits

            plt.xlim(0,1)
            plt.ylim(0,1)
            plt.xlabel("Depth",fontsize=12)
            plt.ylabel("Temperature",fontsize=12)
            plt.tick_params(axis='x',labelsize=10)
            plt.tick_params(axis='y',labelsize=10)

#           Filename and save T_z plot

            full_filename = plot_dir + "T_z_%s.png" %(time_text_file)
            try: os.remove(full_filename)
            except OSError: pass
            plt.savefig(full_filename)
            plt.close()

#-------------------------------------------------------------------------------
#           Plot velocity versus depth
#-------------------------------------------------------------------------------

#           Set the figure size and labels

            fig = plt.figure(figsize=(9,4.5)) ; fig = plt.subplots_adjust(left=0.09,bottom=0.11,top=0.93,right=0.96)

#           Title with time, Ra, Pr, A and benchmark results

            if output_units == "physical": title = r'$V_{rms}$(%s): ' %(time_text_title)
            else                         : title = r'$V_{rms}$(%d): '  %(it)

            Ra_power = np.log10(Ra) ; Ra_int_power = int(Ra_power) ; Ra_rem_power = Ra_power-Ra_int_power ; Ra_factor    = 10**Ra_rem_power
            Pr_power = np.log10(Pr) ; Pr_int_power = int(Pr_power) ; Pr_rem_power = Pr_power-Pr_int_power ; Pr_factor    = 10**Pr_rem_power

            if abs(Ra_factor - 1) < 0.05: title = title + r'$Ra = 10^{%d}$,' %(Ra_int_power)
            else                        : title = title + r'$Ra = %3.1f \times 10^{%d},$' %(Ra_factor,Ra_int_power)

            if abs(Pr_factor - 1) < 0.05: title = title + r'$Pr = 10^{%d},$' %(Pr_int_power)
            else                        : title = title + r'$Pr = %3.1f \times 10^{%d},$' %(Pr_factor,Pr_int_power)
            title = title + r'$A=%5.2f,$' %(A)

            title = title + get_text_from_float(r'$Nu=$',Nu) + "," + get_text_from_float(r'$V_{rms}=$',Vrms) + "," + get_text_from_float(r'$Re=$',Re) + "," + get_text_from_float(r'$Ma=$',Ma)

            plt.title(title,fontsize=12)

            plt.plot   (z,Vz,"-" ,color="blue"                                                  ,linewidth=1.5)

#           Axis labels and plot limits

            plt.xlim(0,1)
            plt.ylim(0,1.02*Vz_max)
            plt.xlabel("Depth",fontsize=12)
            plt.ylabel(r'$V_{rms}$',fontsize=12)
            plt.tick_params(axis='x',labelsize=10)
            plt.tick_params(axis='y',labelsize=10)

#           Filename and save V_z plot

            full_filename = plot_dir + "V_z_%s.png" %(time_text_file)
            try: os.remove(full_filename)
            except OSError: pass
            plt.savefig(full_filename)
            plt.close()

#-------------------------------------------------------------------------------
#           Temperature difference of the TBL relative to T(q_max) versus distance from the boundary and its cubic spline
#-------------------------------------------------------------------------------

            fig = plt.figure(figsize=(7,4.5)) ; fig = plt.subplots_adjust(left=0.12,bottom=0.12,top=0.93,right=0.96)
            title = r'$T$(%s)' %(time_text_title) ; plt.title(title,fontsize=14)

#           Limit plot to the min distance of nmax points and frac x TBL width

            frac   = 2 ; nmax = 8
            numax  = int(TBLW_u_fd+2.5) ; ndmax = int(TBLW_d_fd+2.5)
            dumax  = dx[0]*(numax-1)/L_l ; ddmax = dx[nz-1]*(ndmax-1)/L_l; dmax = max(dumax,ddmax)
            
            xmin   = 0 ; xmax = frac*max(lu_fd,ld_fd) ; xmax = min(xmax,dmax)
            ymin   = 0 ; ymax = 1
            yspace = 0.0 ; yrange = ymax - ymin ; ymin -= yspace*yrange ; ymax += yspace*yrange

            plt.xlim(xmin,xmax) ; plt.ylim(ymin,ymax)

            plt.xlabel(r'Distance from boundary (blue: upper, red: base)',fontsize=12)
            plt.ylabel(r'$T$',fontsize=12)

            plt.plot(z       ,Tz     ,"o-",label=r'$T_u$',color="blue",linewidth=2,markersize=5,markerfacecolor='blue')
            plt.plot(1-z_flip,Tz_flip,"o-",label=r'$T_b$',color="red" ,linewidth=2,markersize=5,markerfacecolor='red' )

            plt.legend(loc="upper left",fontsize=12)

            full_filename = plot_dir + "TBL_T_%s.png" %(time_text_file)
            plt.savefig(full_filename)
            plt.close()

#-------------------------------------------------------------------------------
#           Nu of the TBL versus distance from the boundary and its cubic spline
#-------------------------------------------------------------------------------

            fig = plt.figure(figsize=(7,4.5)) ; fig = plt.subplots_adjust(left=0.12,bottom=0.12,top=0.93,right=0.96)

            title = r'$Nu$(%s) = %5.3g' %(time_text_title,Nu) ; plt.title(title,fontsize=14)

            plt.xlabel(r'Distance from boundary (blue: upper, red: base)',fontsize=12)
            plt.ylabel(r'$Nu$',fontsize=12)

            ymin = np.min(q_fd) ; ymax = np.max(q_fd)
            yspace = 0.05 ; yrange = ymax - ymin ; ymin -= yspace*yrange ; ymax += yspace*yrange
            plt.xlim(xmin,xmax) ; plt.ylim(ymin,ymax)

            plt.plot(z_fd     ,q_fd     ,"-o",label=r'$Nu_u$',color="blue",linewidth=2,markersize=5,markerfacecolor='blue')
            plt.plot(z_fd_flip,q_fd_flip,"-o",label=r'$Nu_d$',color="red" ,linewidth=2,markersize=5,markerfacecolor='red')

            plt.legend(loc="upper right",fontsize=12)

            full_filename = plot_dir + "TBL_Nu_%s.png" %(time_text_file)
            plt.savefig(full_filename)
            plt.close()

#-------------------------------------------------------------------------------
#           Print benchmark summary
#-------------------------------------------------------------------------------

            print("        ======================================",flush=True)
            print("        BENCHMARK SUMMARY: %s" %(time_text_title),flush=True)
            print("        ======================================",flush=True)
            print("        iu        = %5d   = index upper q max" %(iu),flush=True)
            print("        id        = %5d   = index base  q max" %(id),flush=True)
            print("        --------------------------------------",flush=True)
            print("        f = r1/r2 = %9.4g" %(dx[nz-1]/dx[0 ]),flush=True)
            print("        q_cond    = %9.4g = (L/r2)/(ln(r2/r1))" %(q_cond),flush=True)
            print("        --------------------------------------",flush=True)
            print("        lu        = %9.4g = %9.4g lattice cells" %(lu_fd    ,TBLW_u_fd    ),flush=True)
            print("        lb        = %9.4g = %9.4g lattice cells" %(ld_fd    ,TBLW_d_fd    ),flush=True)
            print("        --------------------------------------",flush=True)
            print("        Nu_u      = %9.4g" %(Nu_u_fd)     ,flush=True)
            print("        Nu_b      = %9.4g" %(Nu_d_fd)     ,flush=True)
            print("        ======================================",flush=True)
            print("        Nu        = %9.4g" %(Nu_fd)       ,flush=True)
            print("        Vrms      = %9.4g" %(Vrms)        ,flush=True)
            print("        T_bar     = %9.4g" %(T_bar)       ,flush=True)
            print("        --------------------------------------",flush=True)
            print("        Re        = %9.4g" %(Re)          ,flush=True)
            print("        Ma        = %9.4g" %(Ma)          ,flush=True)
            print("        l         = %9.4g = %9.4g lattice cells" %(l_fd,TBLW_fd),flush=True)
            print("        ======================================",flush=True)

#-------------------------------------------------------------------------------
#       Append benchmark results to benchmark file
#-------------------------------------------------------------------------------

            print("        Appending (bmk) Nu = %8.4g, Vrms = %8.4g Re = %8.4g, Ma = %8.3g" %(Nu,Vrms,Re,Ma),flush=True)
            print("        ======================================",flush=True)

        f_bmk = open(bmk_file,"a")
        f_bmk.writelines("%6d %s %9.4g %9.4g %9.4g %9.4g %9.4g %9.4g %9.4g %9.4g %9.4g %9.4g %9.4g %9.4g %9.4g\n" %(t_bmk,timescale_Inc_bmk,math.log(Ra_l,10),math.log(Pr_l,10),A,Ma,T_bar,Nu_fd,Vrms,Re,TBLW_fd,TBLW_u_fd,TBLW_d_fd,Nu_u_fd,Nu_d_fd))
        f_bmk.close()

