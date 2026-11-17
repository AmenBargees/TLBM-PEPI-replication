#===============================================================================
# Optimal domain decomposition for Lattice Boltzmann code in 2D or 3D allowing
# 1D, 2D or 3D domain decompositions
#
# Author:       Peter Mora
#
# Version date: 20 November 2022
#
# COPYRIGHT
#
# Strictly no distribution allowed outside of the collaboration between
# Peter Mora and Gabriele Morra without permission from the Author, Peter Mora
#===============================================================================

import os
import numpy as np
import time
import matplotlib.pyplot as plt
import matplotlib
matplotlib.rcParams['text.usetex'] = True

#---------------------------------------------------------------------------------------
# Remove Latex stuff from title for print statements
#---------------------------------------------------------------------------------------

def get_print_title(title):
    print_title = ""
    n = 0
    i = 0
    while(i<len(title)):
        if title[i:i+15] == "$\\displaystyle ":
            i += 15
        elif title[i] == "$":
            i += 1
        elif title[i:i+5] == "{\\bf ":
            i += 5
        elif title[i] == "}":
            i += 1
        elif title[i] == "\\":
            i += 1
        else:
            print_title = print_title + title[i]
            n += 1
            i += 1
    print_title = print_title
    return print_title

def plot_field(figure,title,filename,field,mean,clip,gpow,nx,nz):

    global pathnameplots

    mz = field.shape[0]
    mx = field.shape[1]
    tmp = field
    tmp_max = np.max(tmp)
    tmp_min = np.min(tmp)
#    print("Plotting ",title," E [%10.3g" %(tmp_min),",%10.3g" %(tmp_max),"]",flush=True)
    print("Plotting ",get_print_title(title)," E [%10.4g" %(tmp_min),",%10.4g" %(tmp_max),"]",flush=True)

#    print(" abs(tmp_max/tmp_min) = ",np.abs(tmp_max/tmp_min),flush=True)

    tmp_mean = np.sum(tmp/(mx*mz))
#    print("mean = ",tmp_mean,flush=True)

#   If mean = -1, subtract the mean from the field and then set the color scale from -max(abs(field)) to max(abs(field))
#   This makes the plot color scale be the difference between the field and the mean so the mean will always have the same color

    if mean == -1:
        tmp     -= tmp_mean
        tmp_max -= tmp_mean
        tmp_min -= tmp_mean
        tmp_abs_max = max(abs(tmp_max),abs(tmp_min))
        vmin = -tmp_abs_max
        vmax =  tmp_abs_max

#   If mean = 0, assume the mean is zero and set the color scale from -max(abs(field) to max(abs(field))
#   This ensures that a value of zero in the field always has the same color

    elif mean == 0:
        tmp_abs_max = max(abs(tmp_max),abs(tmp_min))
        vmin = -tmp_abs_max
        vmax =  tmp_abs_max

#   Any other value for mean just sets the color scale from min(field) to max(field)
#   This does the usual color scaling from the min(field) to the max(field)

    else:
        vmin = tmp_min
        vmax = tmp_max
    if vmin < 0:
        vmin *= clip
    if vmax > 0:
        vmax *= clip
        if vmin > 0:
            vmin *= clip
    print("vrange = [%10.3g," %(vmin),",%10.3g]" %(vmax),flush=True)

    epsilon = 1.e-30
    abs_tmp = np.abs(tmp)
    max_val = np.max(abs_tmp)
    tmp     = tmp/(max_val+epsilon)
    abs_tmp = abs_tmp/(max_val+epsilon)

    tmp = max_val * tmp * np.power(abs_tmp,gpow)
#    clip_val = clip * max_val
#    tmp = np.clip(tmp,-clip_val,clip_val)

    aspect_field = 1
    aspect_ratio = float(mx)/float(mz)
    if aspect_ratio < 0.75:    # 1x2
        xsize = 3.5
        ysize = 6
#        plt.figure(figure, figsize=(3.5,6)) # 1x2
    elif aspect_ratio < 1.5:   # 1x1
        xsize = 6
        ysize = 5.4
#        plt.figure(figure, figsize=(6,6))   # 1x1
    elif aspect_ratio < 2.5:
        xsize = 6
        ysize = 3.5
        xsize *= 1.1
        ysize *= 1.1
#        plt.figure(figure, figsize=(6,3.5)) # 2x1
    elif aspect_ratio < 3.5:
        xsize = 6
        ysize = 2.5
#        plt.figure(figure, figsize=(6,2.5)) # 3x1
    else:
        xsize = 6
        ysize = 2
#        plt.figure(figure, figsize=(6,2))   # 4x1
#    if high_resolution_plot:
#        xsize *= 2
#        ysize *= 2
    plt.figure(figure, figsize=(xsize,ysize))

    fig = plt.subplots_adjust(left=0.1, bottom=0.05, top=1, right=1.0)
    image = plt.imshow(tmp,extent=[1, mx, mz, 1],aspect=aspect_field,cmap=plt.get_cmap("jet"),vmin=vmin,vmax=vmax)

    plt.title(title,fontsize=14)

    plot_colorbar = True
    if plot_colorbar:
#            norm = mpl.colors.Normalize(vmin=np.min(tmp), vmax=np.max(tmp))
            where = "Saudi"
            if aspect_ratio <= 4:
                if where == "Saudi":
                    cbar = plt.colorbar(image,orientation='vertical', shrink=0.68, pad=0.05,aspect=24/aspect_ratio)
                else:
                    cbar = plt.colorbar(image,orientation='vertical', shrink=0.82, pad=0.05, aspect=14/aspect_ratio)
            else:
                cbar = plt.colorbar(image,orientation='vertical', shrink=1.3, pad=0.05, aspect=2*14/aspect_ratio)

    if figure == 1:
        plt.xlim(1, mx)
        plt.ylim(mz, 1)
    else:
        print("xlim = 1:65, ylim = ",nz,":",nz-65)
        plt.xlim(1, 65)
        plt.ylim(65, 1)
    plt.tick_params(axis='x',labelsize=12)
    plt.tick_params(axis='y',labelsize=12)
    plt.xlabel("nx_d",fontsize=16)
    plt.ylabel("nz_d",fontsize=16)
    full_filename = pathnameplots + filename + '.png'
    try: os.remove(full_filename)
    except OSError: pass
    print("filename = ",full_filename,flush=True)
    plt.savefig(full_filename)
    plt.close()


#=======================================================================
# Domain decomposition in 1D, 2D or 3D
#=======================================================================

def my_domain_decomposition(size,rank,nx,ny,nz,nxi_min,nxi_max,nyi_min,nyi_max,nzi_min,nzi_max,domain_decomposition):



#   If verbose = True, get verbose prints showing progress

    global verbose

    verbose = False

    if   domain_decomposition == "1D":
        nyi_min = ny ; nyi_max = ny
        nzi_min = nz ; nzi_max = nz
        nr,nrx,nry,nrz,nxi,nyi,nzi,nx_i,ny_i,nz_i,mx,my,mz = domain_decomposition_ND(size,rank,nx,1,nz,nxi_min,nxi_max,nyi_min,nyi_max,nzi_min,nzi_max)
    elif domain_decomposition == "2D":
        nyi_min = ny ; nyi_max = ny
        nr,nrx,nry,nrz,nxi,nyi,nzi,nx_i,ny_i,nz_i,mx,my,mz = domain_decomposition_ND(size,rank,nx,1,nz,nxi_min,nxi_max,nyi_min,nyi_max,nzi_min,nzi_max)
    elif domain_decomposition == "3D":
        nr,nrx,nry,nrz,nxi,nyi,nzi,nx_i,ny_i,nz_i,mx,my,mz = domain_decomposition_ND(size,rank,nx,ny,nz,nxi_min,nxi_max,nyi_min,nyi_max,nzi_min,nzi_max)

    if rank == 0 and verbose:
        if nrx*nry*nrz != size:
            print("NOTE: For 100% core useage specify ",nrx*nry*nrz," cores",flush=True)
            print("============================================",flush=True)

    return nr,nrx,nry,nrz,nxi,nyi,nzi,nx_i,ny_i,nz_i,mx,my,mz

#-------------------------------------------------------------------------------
# Get the 2D or 3D domain decomposition
#-------------------------------------------------------------------------------

def domain_decomposition_ND(size,rank,nx,ny,nz,nxi_min,nxi_max,nyi_min,nyi_max,nzi_min,nzi_max):

    global coreUsed

#   Get the number of domains and domain sizes

    (nrx,nry,nrz,nxi,nyi,nzi,nxlast,nylast,nzlast)=core_fit_ND(size,rank,nx,ny,nz,nxi_min,nxi_max,nyi_min,nyi_max,nzi_min,nzi_max)

#   Set number of cores used and max rank

    coreUsed=nrx*nry*nrz
    nr = nrx*nry*nrz

#   Add space for edges if needed

    if nrx == 1: nx_i = nxi
    else:        nx_i = nxi+2
    if nry == 1: ny_i = nyi
    else:        ny_i = nyi+2
    if nrz == 1: nz_i = nzi
    else:        nz_i = nzi+2

#   Set x, y and z domain sizes

    mx = np.ones(nrx,dtype=int)*nxi
    my = np.ones(nry,dtype=int)*nyi
    mz = np.ones(nrz,dtype=int)*nzi
    if nrx > 1: mx[-1]=nxlast
    if nry > 1: my[-1]=nylast
    if nrz > 1: mz[-1]=nzlast

    return nr,nrx,nry,nrz,nxi,nyi,nzi,nx_i,ny_i,nz_i,mx,my,mz

#===============================================================================
# Get the optimal ND domain decomposition that best fits onto cores
#===============================================================================

#-------------------------------------------------------------------------------
# Function to compute the edge domain sizes
#-------------------------------------------------------------------------------

def edge_domain_sizes(nx,ny,nz,nd_x,nd_y,nd_z,mx,my,mz):

    nx_sum  = (nd_x-1)*mx
    nx_last =  nx - nx_sum
    ny_sum  = (nd_y-1)*my
    ny_last =  ny - ny_sum
    nz_sum  = (nd_z-1)*mz
    nz_last =  nz - nz_sum

    return  nx_last,ny_last,nz_last

#-------------------------------------------------------------------------------
# Function to compute the calculation cost of a domain decomposition
#-------------------------------------------------------------------------------

def calculation_cost(nd_x,nd_y,nd_z,mx,my,mz,nx,ny,nz):

    x_calc_cost    = nd_x*mx
    y_calc_cost    = nd_y*my
    z_calc_cost    = nd_z*mz

    calc_cost      = (x_calc_cost * y_calc_cost * z_calc_cost)/(nx*ny*nz)

    return calc_cost

#-------------------------------------------------------------------------------
# Function to calculate communications cost of overheads
# and for data transfer for a domain decomposition
#-------------------------------------------------------------------------------

def communication_cost(nd_x,nd_y,nd_z,mx,my,mz):

#   Communications overheads proportional to number of sends

    if nd_x <= 1.001: x_comm = 0
    else:         x_comm = 1
    if nd_y <= 1.001: y_comm = 0
    else:         y_comm = 1
    if nd_z <= 1.001: z_comm = 0
    else:         z_comm = 1

    comm_oh_cost       = nd_x*x_comm + nd_y*y_comm + nd_z*z_comm

#   Communications data transfer proportional to volume of data sent allowing
#   for the need to send extra data if there are edges to domains

    kx = 2
    ky = 2
    kz = 2

    if nd_x <= 1.001: kx = 0
    if nd_y <= 1.001: ky = 0
    if nd_z <= 1.001: kz = 0

    mx_d = mx+kx
    my_d = my+ky
    mz_d = mz+kz

    comm_data_cost       = my_d*mz_d*x_comm + mx_d*mz_d*y_comm + mx_d*my_d*z_comm

    return comm_oh_cost,comm_data_cost

def core_fit_ND(size,rank,nx,ny,nz,nxi_min,nxi_max,nyi_min,nyi_max,nzi_min,nzi_max):

    global pathnameplots

#   Set the communications cost parameters to be used int he optimization

    rel_comm    = 0.5   # Communications cost relative to calculations (at linear speedup breakdown)
    rel_comm_oh = 0.05  # Assumed communications overhead (comm overhead relative to the comms)

#   Set the theoretical optimal domain decomposition => square or cubic and exactly fitting onto mesh

    if ny == 1:
        ix    = np.sqrt(nx*nz/size)
        iy    = ny
        iz    = ix
        id_x  = nx/ix
        id_y  = 1
        id_z  = nz/iz
#        if rank == 0:
#            print("nx = ",nx," <=? id_x*ix = ",id_x*ix,", id_x = ",id_x,", ix = ",ix,flush=True)
#            print("nz = ",nz," <=? id_z*iz = ",id_z*iz,", id_z = ",id_z,", iz = ",iz,flush=True)
    else:
        if nxi_min == nxi_max and nxi_min == nx:
            ix    = nx
            iy    = np.sqrt(ny*nz/size)
            iz    = iy
            id_x  = 1
            id_y  = ny/iy
            id_z  = nz/iz
#            if rank == 0:
#                print("ny = ",ny," <=? id_y*iy = ",id_y*iy,", id_y = ",id_y,", iy = ",iy,flush=True)
#                print("nz = ",nz," <=? id_z*iz = ",id_z*iz,", id_z = ",id_z,", iz = ",iz,flush=True)
        else:
            id_x  = np.power(size,0.33333)
            id_y  = id_x ; id_z = id_x
            ix    = nx/id_x
            iy    = ny/id_y
            iz    = nz/id_z

            ix    = np.power(nx*ny*nz/size,0.33333)
            iy    = ix
            iz    = ix
            id_x  = nx/ix
            id_y  = ny/iy
            id_z  = nz/iz
#            if rank == 0:
#                print("nx = ",nx," <=? id_x*ix = ",id_x*ix,", id_x = ",id_x,", ix = ",ix,flush=True)
#                print("ny = ",ny," <=? id_y*iy = ",id_y*iy,", id_y = ",id_y,", iy = ",iy,flush=True)
#                print("nz = ",nz," <=? id_z*iz = ",id_z*iz,", id_z = ",id_z,", iz = ",iz,flush=True)

    nd_x_optimal = id_x
    nd_y_optimal = id_y
    nd_z_optimal = id_z
    mx_optimal   = ix
    my_optimal   = iy
    mz_optimal   = iz

    calc_cost_optimal                           = calculation_cost(id_x,id_y,id_z,ix,iy,iz,nx,ny,nz)
    comm_cost_oh_optimal,comm_cost_data_optimal = communication_cost(id_x,id_y,id_z,ix,iy,iz)

#   Set ranges to scan for optimal region size per node but ensure valid ranges

#   Check z scan range

    nxi1 = min(nxi_min,nxi_max) ; nxi2 = max(nxi_min,nxi_max)
    nxi_min = nxi1 ; nxi_max = nxi2
    
    m_min_x = min(nx,nxi_min) ; m_max_x = min(nx,nxi_max)

    if size*m_min_x < nx: m_min_x = max(1,nx//size)
    if m_max_x < m_min_x: m_max_x = m_min_x

#   Check y scan range

    if nyi_min == "unspecified": nyi_min = 1 
    if nyi_max == "unspecified": nyi_max = ny

    nyi1 = min(nyi_min,nyi_max) ; nyi2 = max(nyi_min,nyi_max)
    nyi_min = nyi1 ; nyi_max = nyi2

    m_min_y = min(ny,nyi_min) ; m_max_y = min(ny,nyi_max)

    if size*m_min_y < ny: m_min_y = max(1,ny//size)
    if m_max_y < m_min_y: m_max_y = m_min_y

#   Check z scan range

    if nzi_min == "unspecified": nzi_min = 1 
    if nzi_max == "unspecified": nzi_max = nz

    nzi1 = min(nzi_min,nzi_max) ; nzi2 = max(nzi_min,nzi_max)
    nzi_min = nzi1 ; nzi_max = nzi2

    m_min_z = min(nz,nzi_min) ; m_max_z = min(nz,nzi_max)

    if size*m_min_z < nz: m_min_z = max(1,nz//size)
    if m_max_z < m_min_z: m_max_z = m_min_z

#   Loop over possible domain sizes to obtain the optimal domain size per node

    if rank == 0 and verbose:
        print("===================================================================",flush=True)
        print("Domain decomposition to optimize efficiency on ",size," cores",flush=True)
        print("for a model size of (",nx,",",ny,",",nz,")",flush=True)
        print("-------------------------------------------------------------------",flush=True)
        print("Optimization takes into account core usage, computational cost",flush=True)
        print("and communication overhead and the time to send data",flush=True)
        print("===================================================================",flush=True)
        print("Assumptions:",flush=True)
        print("------------",flush=True)
        print("   (1) (comm cost)/(calc cost) = ",100*rel_comm," %, and",flush=True)
        print("   (2) (comms oh cost)/(comm send cost) = ",100*rel_comm_oh," %",flush=True)
        print()
        print("Note: efficiencies are relative to the",flush=True)
        print("theoretical optimal domain decomposition",flush=True)
        print("-------------------------------------------------------------------",flush=True)
        print("THEORETICAL OPTIMAL DOMAIN DECOMPOSITION",flush=True)
        print("-------------------------------------------------------------------",flush=True)
        print("                     x        y        z",flush=True)
        print("Model size:   = %6d   %6d   %6d" %(nx,ny,nz),flush=True)
        print("-------------------------------------------------------------------",flush=True)
        print("# domains     = %6.1f   %6.1f   %6.1f" %(nd_x_optimal,nd_y_optimal,nd_z_optimal),flush=True)
        print("Domain sizes  = %6.1f   %6.1f   %6.1f" %(mx_optimal,my_optimal,mz_optimal),flush=True)
        print("===================================================================",flush=True)
        print("STARTING DOMAIN DECOMPOSITION OPTIMIZATION LOOPS",flush=True)
        print("-------------------------------------------------------------------",flush=True)
        print("Domain size scan x range is set from ",m_min_x," to ",m_max_x,flush=True)
        print("Domain size scan y range is set from ",m_min_y," to ",m_max_y,flush=True)
        print("Domain size scan z range is set from ",m_min_z," to ",m_max_z,flush=True)
        print("-------------------------------------------------------------------",flush=True)

#---------------------------------------------------------
#   Start optimization loop
#---------------------------------------------------------

    total_efficiency  = 0
    cores_used        = 0
    ix_print_inc = max(1,(m_max_x-m_min_x)//10)
    iy_print_inc = max(1,(m_max_y-m_min_y)//10)
    iz_print_inc = max(1,(m_max_z-m_min_z)//10)

    last_print = False

#   Option to be able to plot the objective function for DD (only for 2D DD)

    plot_DD = False  # For 2D DD can set to True to get the plot of DD efficiency vs x and z

    if plot_DD == True and ny == 1: efficiency = np.zeros((nz,nx))

    for ix in range(m_min_x,m_max_x+1):

#      Get number of x domains

       id_x = int(nx//ix)
       if id_x*ix < nx: id_x += 1

       if verbose:
           if rank == 0 and (ix == m_min_x or (ny == 1 and ix%ix_print_inc == 0)):
               if last_print: print("--------------------------------------------",flush=True)
               print("x domain scan %8d      to %8d" %(ix,min(m_max_x,min(ix+ix_print_inc,nx))),flush=True)
           last_print = False


#      Save efficiency as a function of the y and z domain sizes to be able to break from z loop

       efficiency_yz = np.zeros((m_max_y,m_max_z))

       for iy in range(m_min_y,m_max_y+1):

          if rank == 0 and not verbose and ny > 1 and ix == m_min_x and iy%iy_print_inc == 0: print(".",end="",flush=True)

#         Get number of y domains

          id_y = int(ny//iy)
          if id_y*iy < ny: id_y += 1

          for iz in range(m_min_z,m_max_z+1):

                if rank == 0 and not verbose and ny == 1 and ix == m_min_x and iz%iz_print_inc == 0: print(".",end="",flush=True)


#               Get number of z domains

                id_z = int(nz//iz)
                if id_z*iz < nz: id_z += 1

#               Do optimization work only for viable domain sizes

                icores_used         = id_x*id_y*id_z
                if icores_used <= size:

#                   Calculate the usual and last x and z size of the domains

                    ix_last,iy_last,iz_last = edge_domain_sizes(nx,ny,nz,id_x,id_y,id_z,ix,iy,iz)

#                   Core efficiency

                    icore_efficiency      = icores_used/size

#                   Computational efficiency

                    icalc_cost       = calculation_cost(id_x,id_y,id_z,ix,iy,iz,nx,ny,nz)
                    icalc_efficiency = calc_cost_optimal/icalc_cost

#                   Communications efficiency allowing for cost of overhead versus cost of data transfer

                    icomm_oh_cost,icomm_data_cost = communication_cost(id_x,id_y,id_z,ix,iy,iz)

                    if icomm_oh_cost == 0:
                        icomm_oh_efficiency = 1
                        icomm_data_efficiency = 1
                    else:
                        icomm_oh_efficiency   = comm_cost_oh_optimal  /icomm_oh_cost
                        icomm_data_efficiency = comm_cost_data_optimal/icomm_data_cost

                    icomm_efficiency      = rel_comm_oh*icomm_oh_efficiency +(1-rel_comm_oh)*icomm_data_efficiency

#                   Total efficiency

                    itotal_efficiency  = icore_efficiency*((1-rel_comm)*icalc_efficiency + rel_comm*icomm_efficiency)

#                   Save data for plot if needed

                    if plot_DD: efficiency[iz-1,ix-1] = itotal_efficiency

#                   Optimize to maximize core utilization efficiency and reduce the total cost per core

                    if (itotal_efficiency > total_efficiency):
                        improvement_efficiency = itotal_efficiency-total_efficiency
                        improvement_cores_used = icores_used - cores_used
                        cores_used             = icores_used
                        core_efficiency        = icore_efficiency
                        calc_efficiency        = icalc_efficiency
                        calc_cost              = icalc_cost
                        comm_efficiency        = icomm_efficiency
                        comm_oh_efficiency     = icomm_oh_efficiency
                        comm_data_efficiency   = icomm_data_efficiency
                        total_efficiency       = itotal_efficiency
                        nd_x    = id_x
                        nd_y    = id_y
                        nd_z    = id_z
                        mx      = ix
                        my      = iy
                        mz      = iz
                        mx_last = ix_last
                        my_last = iy_last
                        mz_last = iz_last
                        last_print = True
                        if rank == 0 and verbose:
                            print("============================================",flush=True)
                            print("Accepting new domain decomposition",flush=True)
                            print("============================================",flush=True)
                            print("                    x       y       z",flush=True)
                            print("Model size:   = %5d   %5d   %5d" %(nx,ny,nz),flush=True)
                            print("--------------------------------------------",flush=True)
                            print("# domains     = %5d   %5d   %5d" %(nd_x,nd_y,nd_z),flush=True)
                            print("Domain sizes  = %5d   %5d   %5d" %(mx,my,mz),flush=True)
                            if mx_last != mx or my_last != my or mz_last != mz:
                                print("Edge   sizes = %5d , %5d , %5d" %(mx_last,my_last,mz_last),flush=True)
                                print("--------------------------------------------",flush=True)
                                print("Total size    = %5d   %5d   %5d" %((nd_x-1)*mx+mx_last,(nd_y-1)*my+my_last,(nd_z-1)*mz+mz_last),flush=True)
                                print("--------------------------------------------",flush=True)
                                print("Edge domains not full, sub-optimal efficiency",flush=True)
                            else:
                                print("Total size    = %5d   %5d   %5d" %(nd_x*mx,nd_y*my,nd_z*mz),flush=True)
                                print("--------------------------------------------",flush=True)
                                print("All domains full, optimal efficiency",flush=True)
                            print("--------------------------------------------",flush=True)
                            print("Total cores              = ",size,flush=True)
                            print("Used  cores              = ",nd_x*nd_y*nd_z,flush=True)
                            print("Idle  cores              = ",size-nd_x*nd_y*nd_z,flush=True)
                            print("--------------------------------------------",flush=True)
                            print("Core usage               = %6.2f" %(100*cores_used/size)," %",flush=True)
                            print("--------------------------------------------",flush=True)
                            print("Calculation   efficiency = %6.2f" %(100*calc_efficiency)," %",flush=True)
                            print("Communication efficiency = %6.2f" %(100*comm_efficiency)," %",flush=True)
                            print("--------------------------------------------",flush=True)
                            print("comm oh efficiency       = %6.2f" %(100*comm_oh_efficiency)," %",flush=True)
                            print("comm data efficiency     = %6.2f" %(100*comm_data_efficiency)," %",flush=True)
                            print("--------------------------------------------",flush=True)
                            print("Total efficiency         = %6.2f" %(100*total_efficiency)," %",flush=True)
                            print("--------------------------------------------",flush=True)
                            print("Core usage improvement   = %7.3f" %(100*improvement_cores_used/size)," %",flush=True)
                            print("Efficiency improvement   = %7.3f" %(100*improvement_efficiency)," %",flush=True)

#                   Exit z scan after a max is achieved

                    efficiency_yz[iy-1,iz-1] = itotal_efficiency
                    if not plot_DD and \
                           (ny == 1 and efficiency_yz[iy-1,iz-1] < efficiency_yz[iy-1,iz-2]) or \
                           (ny >  1 and efficiency_yz[iy-1,iz-1] < efficiency_yz[iy-2,iz-1]  \
                                    and efficiency_yz[iy-1,iz-1] < efficiency_yz[iy-1,iz-2]):
                        break

    if plot_DD and ny == 1:
        pathnameplots = "./"
        plot_field(1,"DD efficiency for a "+str(nx)+"x"+str(nz)+" lattice on "+str(size)+" cores","Efficiency_DD_"+str(nx)+"x"+str(nz)+"_"+str(size),100*efficiency[0:min(nz,65),:min(nx,65)],1,1,1,nx,nz)

    if rank == 0 and verbose:
        print("============================================",flush=True)
        print("OPTIMAL DOMAIN DECOMPOSITION",flush=True)
        print("============================================",flush=True)
        print("                    x       y       z",flush=True)
        print("Model size:   = %5d   %5d   %5d" %(nx,ny,nz),flush=True)
        print("--------------------------------------------",flush=True)
        print("# domains     = %5d   %5d   %5d" %(nd_x,nd_y,nd_z),flush=True)
        print("Domain sizes  = %5d   %5d   %5d" %(mx,my,mz),flush=True)
        if mx_last != mx or my_last != my or mz_last != mz:
            print("Edge    sizes = %5d   %5d   %5d" %(mx_last,my_last,mz_last),flush=True)
            print("--------------------------------------------",flush=True)
            print("Total size    = %5d   %5d   %5d" %((nd_x-1)*mx+mx_last,(nd_y-1)*my+my_last,(nd_z-1)*mz+mz_last),flush=True)
            print("--------------------------------------------",flush=True)
            print("Edge domains not full, sub-optimal efficiency",flush=True)
        else:
            print("Total size    = %5d   %5d   %5d" %(nd_x*mx,nd_y*my,nd_z*mz),flush=True)
            print("--------------------------------------------",flush=True)
            print("All domains full, optimal efficiency",flush=True)

        print("--------------------------------------------",flush=True)
        print("Total cores              = ",size,flush=True)
        print("Used  cores              = ",nd_x*nd_y*nd_z,flush=True)
        print("Idle  cores              = ",size-nd_x*nd_y*nd_z,flush=True)
        print("--------------------------------------------",flush=True)
        print("Core usage               = %6.2f" %(100*cores_used/size)," %",flush=True)
        print("--------------------------------------------",flush=True)
        print("Calculation   efficiency = %6.2f" %(100*calc_efficiency)," %",flush=True)
        print("--------------------------------------------",flush=True)
        print("Communication efficiency = %6.2f" %(100*comm_efficiency)," %",flush=True)
        print("--------------------------------------------",flush=True)
        print("comm oh efficiency       = %6.2f" %(100*comm_oh_efficiency)," %",flush=True)
        print("comm data efficiency     = %6.2f" %(100*comm_data_efficiency)," %",flush=True)
        print("--------------------------------------------",flush=True)
        print("Total efficiency         = %6.2f" %(100*total_efficiency)," %",flush=True)
        if (mx_last != mx or my_last != my or mz_last != mz):
            if verbose:
                print("============================================",flush=True)
                print("FOR OPTIMAL PERFORMANCE on ",nd_x*nd_y*nd_z," CORES",flush=True)
                print("with full edge domains, use a model size of",flush=True)
                print("============================================",flush=True)
                print("                    x       y       z",flush=True)
                print("Model size:   = %5d   %5d   %5d" %(nd_x*mx,nd_y*my,nd_z*mz),flush=True)
                print("--------------------------------------------",flush=True)
                print("# domains     = %5d   %5d   %5d" %(nd_x,nd_y,nd_z),flush=True)
                print("Domain sizes  = %5d   %5d   %5d" %(mx,my,mz),flush=True)
                print("--------------------------------------------",flush=True)
                print("Total size    = %5d   %5d   %5d" %(nd_x*mx,nd_y*my,nd_z*mz),flush=True)
                print("--------------------------------------------",flush=True)
                print("All domains full, optimal efficiency",flush=True)
                print("--------------------------------------------",flush=True)
            print("============================================",flush=True)

#           Core efficiency

            icore_efficiency = nd_x*nd_y*nd_z/size

#           Computational efficiency (100% by definition)

            icalc_cost       = calculation_cost(nd_x,nd_y,nd_z,mx,my,mz,nx,ny,nz)
            icalc_efficiency = 1

#           Communications efficiency allowing for cost of overhead versus cost of data transfer

            icomm_oh_cost,icomm_data_cost = communication_cost(nd_x,nd_y,nd_z,mx,my,mz)
            if icomm_oh_cost == 0:
                icomm_oh_efficiency = 1
            else:
                icomm_oh_efficiency   = comm_cost_oh_optimal  /icomm_oh_cost
            icomm_data_efficiency = comm_cost_data_optimal/icomm_data_cost

            icomm_efficiency      = rel_comm_oh*icomm_oh_efficiency +(1-rel_comm_oh)*icomm_data_efficiency

#           Total efficiency

            itotal_efficiency  = icore_efficiency*((1-rel_comm)*icalc_efficiency + rel_comm*icomm_efficiency)

            if verbose:
                print("Core usage               = %6.2f" %(100*cores_used/size)," %",flush=True)
                print("--------------------------------------------",flush=True)
                print("Calculation   efficiency = %6.2f" %(100*icalc_efficiency)," %",flush=True)
                print("Communication efficiency = %6.2f" %(100*icomm_efficiency)," %",flush=True)
                print("--------------------------------------------",flush=True)
                print("comm oh efficiency       = %6.2f" %(100*icomm_oh_efficiency)," %",flush=True)
                print("comm data efficiency     = %6.2f" %(100*icomm_data_efficiency)," %",flush=True)
                print("--------------------------------------------",flush=True)
                print("Total efficiency         = %6.2f" %(100*itotal_efficiency)," %",flush=True)
                print("--------------------------------------------")

    return (nd_x,nd_y,nd_z,mx,my,mz,mx_last,my_last,mz_last)

