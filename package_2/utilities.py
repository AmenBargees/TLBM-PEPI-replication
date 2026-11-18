import numpy as np
import sys
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

# Initialize conversion factors

sec_m    = 60              # Seconds to Mins        conversion factor
sec_h    = 60     * sec_m  # Seconds to Hours       conversion factor
sec_d    = 24     * sec_h  # Seconds to Days        conversion factor
sec_w    = 7      * sec_d  # Seconds to Weeks       conversion factor

sec_y    = 365.25 * sec_d  # Seconds to Years       conversion factor
sec_M    = sec_y/12        # Seconds to Months      conversion factor

sec_dy   = sec_y  * 10     # Seconds to Deca  Years conversion factor
sec_hy   = sec_y  * 100    # Seconds to Hecto Years conversion factor
sec_ky   = sec_y  * 1000   # Seconds to Ky          conversion factor

sec_dky  = sec_ky * 10     # Seconds to Deca  Ky    conversion factor
sec_hky  = sec_ky * 100    # Seconds to Hecto Ky    conversion factor
sec_my   = sec_ky * 1000   # Seconds to My          conversion factor

sec_dmy  = sec_my * 10     # Seconds to Deca  My    conversion factor
sec_hmy  = sec_my * 100    # Seconds to Hecto My    conversion factor
sec_gy   = sec_my * 1000   # Seconds to Gy          conversion factor

sec_dgy  = sec_gy * 10     # Seconds to Deca  Gy    conversion factor
sec_hgy  = sec_gy * 100    # Seconds to Hecto Gy    conversion factor
sec_ty   = sec_gy * 1000   # Seconds to Ty          conversion factor

u_cm_py  = 100   *sec_y    # m/sec to cm/y          conversion factor

#-------------------------------------------------------------------------------
# Timescales and conversion factors
#-------------------------------------------------------------------------------

valid_timescales    = ["timesteps", "Ktimesteps" , "Mtimesteps", "Secs", "Mins" ,"Hours", "Days" , "Weeks", "Months", "Years", "dY"   , "hY"     , "Ky"  , "dKy"   , "hKy"  , "My"  , "dMy"  , "hMy"  , "Gy"  , "dGy"  , "hGy"  , "Ty"  ]
factor_timesteps    = [ 1         ,  1000        ,  1000000    ,  1     , sec_m , sec_h ,  sec_d , sec_w  ,  sec_M  ,  sec_y ,  sec_dy, sec_hy   , sec_ky, sec_dky , sec_hky, sec_my, sec_dmy, sec_hmy, sec_gy, sec_dgy, sec_hgy, sec_ty]
physical_timescales = valid_timescales                           [3:]
conversion_timesteps_to_physical_time = factor_timesteps         [3:]

#-------------------------------------------------------------------------------
# Get a reasonable int timevalue and timescale from a float time increment
#-------------------------------------------------------------------------------

def get_int_timevalue_and_timescale_from_time(time,dt_p,rank=0):

    reasonable_timescales = ["Gy","My","Ky","Years","Months","Weeks","Days",",Hours","Secs","Millisecs"]

    found_timescale = False
    for i in range(len(reasonable_timescales)):
        for j in range(len(valid_timescales)):
            time_float     = (time/dt_p)*factor_timesteps[i]
            time_int       = int(time_float/1000+0.5)
            if reasonable_timescales[i] == valid_timescales[j] and time_int >= 1:
                timescale      = valid_timescales[j]
                if rank == 0:
                    print("reasonable_timescales[%d] = %s == %s = valid_timescales[%d]" %(i,reasonable_timescales[i],valid_timescales[j],j),file=sys.stderr,flush=True)
                    print("SELECTED: time_float = %8.3g, time_int = %d, timescale = %s" %(time_float,time_int,timescale),file=sys.stderr,flush=True)
                found_timescale = True
                break
        if found_timescale: break

    return time_int, timescale

#-------------------------------------------------------------------------------
# Convert an input time value in time scale units to timesteps
#-------------------------------------------------------------------------------

def timevalue_to_timesteps(value,timescale,dt_p,rank=0):

    check_timevalue_to_timesteps = False # Option to test time_value_to_timesteps is correctly calculated

    for i in range(len(valid_timescales)):
        if timescale == valid_timescales[i]:
            factor    = factor_timesteps[i]
            if check_timevalue_to_timesteps:
                if rank == 0:
                    for iout in range(2):
                        if   iout == 0: file_out = sys.stderr
                        elif iout == 1: file_out = sys.stdout
                        print("--------------------------------------------------",file=file_out,flush=True)
                        if value == "default" or value == "unspecified" or value == "calculate":
                            print("value = %s" %(value)                           ,file=file_out,flush=True)
                        else:
                            print("value = %s = "    %(value   )           ,end="",file=file_out,flush=True)
                            if   isinstance(value,int  ): print("int"      ,end="",file=file_out,flush=True)
                            elif isinstance(value,float): print("float"    ,end="",file=file_out,flush=True)
                            else                        : print("invalid"  ,end="",file=file_out,flush=True)
                            print(", timescale = %s" %(timescale)                 ,file=file_out,flush=True)
            timesteps_float = value*factor/dt_p
            timesteps_int   = int(timesteps_float+0.5)
            break
        else:
            timesteps_float = 0
            timesteps_int   = 0

    return timesteps_int, timesteps_float, factor

#----------------------------------------------------------------------------
# Get the time text for output files
#----------------------------------------------------------------------------

def timetext(timesteps,timescale,dt_p):

    for i in range(len(valid_timescales)):
        if   timescale == "timesteps":
            text = "%010d" %(timesteps)
        elif timescale == valid_timescales[i]:
            it   = int(dt_p*timesteps/factor_timesteps[i]+0.5)
            text = "%010d" %(it) + "_%s" %(timescale)

    return text

#----------------------------------------------------------------------------
# Get the time value and timescale from the input timetext
#----------------------------------------------------------------------------

def timevalue(parameter,timetext,error_timescale=False,rank=0):
 
    parts = timetext.split("_")

    timevalue = parts[0]
    if len(parts) > 1:
        timescale = parts[1]
    else:
        timescale = "unknown"

    error_timescale = False

    check_timescale = False
    if check_timescale:
        for iout in range(2):
            if   iout == 0: file_out = sys.stderr
            elif iout == 1: file_out = sys.stdout
            print("parameter = %s, timetext = %s, " %(parameter,timetext),end="",file=file_out,flush=True)
            print("parts = ",parts                                              ,file=file_out,flush=True)

    error_timescale = True
    for i in range(len(valid_timescales)):
        if timevalue == "default" or timevalue == "calculate" or timevalue == "unspecified":
            timevalue = 0
            error_timescale = False
            break
        elif len(parts) == 2 and timescale == valid_timescales[i]:
            timevalue = int(timevalue)
            timescale = valid_timescales[i]
            error_timescale = False
            break

    if error_timescale and rank == 0:
        for iout in range(2):
            if   iout == 0: file_out = sys.stderr
            elif iout == 1: file_out = sys.stdout
            print("ERROR: Invalid timescale \"%s\" in time specification \"%s = %s\". Time specifications must be of form \"n_timescale\"" %(timescale,parameter,timetext),file=file_out,flush=True)
            print("       where n is an integer number in units of timescale which must be one of:"       ,file=file_out,flush=True)
            print("           "                                                                                ,end="",file=file_out,flush=True)
            for j in range(len(valid_timescales)):
                if j > 0: print(", ",end="",file=file_out,flush=True)
                print("%s" %(valid_timescales[j]),end="",file=file_out,flush=True)
            print(file=file_out,flush=True)

            if check_timescale and error_timescale:
                print("       timevalue       = ",timevalue,file=file_out,flush=True)
                print("       timescale       = ",timescale,file=file_out,flush=True)
                print("       error_timescale = ",error_timescale,file=file_out,flush=True)

    return timevalue, timescale, error_timescale

#----------------------------------------------------------------------------
# Get seconds from time text in format n_timescale
#----------------------------------------------------------------------------

def timevalue_to_secs(timevalue,timescale,rank=0):

    if timescale != "unspecified":
        for i in range(len(valid_timescales)):
            if timescale == valid_timescales[i]:
                secs = factor_timesteps[i]*timevalue
    else:
        secs = float(timevalue)

    return secs

#----------------------------------------------------------------------------
#   Get the time text to print from a time specification
#----------------------------------------------------------------------------

def get_time_value_text_int_timescale_from_time_spec(time_spec,time_value,dt_p):

    error = False
    time_int,time_scale,error    = timevalue(time_spec,time_value,error,0)
    time_steps,time_steps_float,time_steps_factor = timevalue_to_timesteps(time_int,time_scale,dt_p,0)

    return time_steps, time_int, time_scale

#----------------------------------------------------------------------------
# Get time specification and its integer value & timescale from a time in secs
#----------------------------------------------------------------------------

def time_to_time_spec(time_secs,rank):

    check_function = False

    for i in range(len(physical_timescales)-1,-1,-1):
        t_p = time_secs/conversion_timesteps_to_physical_time[i]
        if rank == 0 and check_function:
            print("physical_timescales[%d] = %s, t_p = %8.3g" %(i,physical_timescales[i],t_p),flush=True)
        if t_p >= 1:
            time_float     = t_p
            time_int       = int(time_float)
            timescale_time = "%s" %(physical_timescales[i])
            time_spec = "%d_%s" %(time_int,timescale_time)
            if rank == 0 and check_function:
                print("t_p >= 0, time_spec = %s, time_int = %d, time_float = %4.2f, timescale_time = %s" %(time_spec,time_int,time_float,timescale_time),flush=True)
    return time_spec, time_int, time_float, timescale_time

#----------------------------------------------------------------------------
# Get integer time from time in secs for a specified timescale
#----------------------------------------------------------------------------

def secs_to_float_time_timescale(secs,time_scale,rank):

    for i in range(len(physical_timescales)):
        if time_scale == physical_timescales[i]:
            t_p = secs/conversion_timesteps_to_physical_time[i]
    return t_p

#-------------------------------------------------------------------------------
# Get the time from the timestep for prints and plots
#
# Arguments are t (LBM time) and dt_p (physical time step) and output_units (physical or timesteps)
#-------------------------------------------------------------------------------

def timesteps_to_time_text(t,dt_p,output_units="physical"):

    t_p = t*dt_p    # Get time in secs

    if output_units == "physical": time_text = time_to_physical_time_text(t_p,dt_p)
    else                         : time_text = "%d timesteps" %(int(t+0.5))

    return time_text

def time_to_physical_time_text(t_p,dt_p):

    if t_p < sec_h:       # Less than 1 hour -> secs
        physical_time_text = "%5.3f" %(t_p)        + " secs"
    elif t_p < sec_d:     # Less than 1 day -> hours
        physical_time_text = "%5.3f" %(t_p/sec_h)  + " hours"
    elif t_p < sec_y:     # Less than 1 year -> days
        physical_time_text = "%5.3f" %(t_p/sec_d)  + " days"
    elif t_p < sec_ky:    # Less than 1 ky   -> years
        physical_time_text = "%5.3f" %(t_p/sec_y)  + " years"
    elif t_p < sec_my:    # Less than 1 my   -> ky
        physical_time_text = "%5.3f" %(t_p/sec_ky) + " Ky"
    elif t_p < sec_gy:    # Less than 1 gy   -> my
        physical_time_text = "%5.3f" %(t_p/sec_my) + " My"
    else:                 # >=        1 my   -> gy
        physical_time_text = "%5.3f" %(t_p/sec_gy) + " Gy"

    return physical_time_text

#-------------------------------------------------------------------------------
# Get a specified plane and coordinate of a 3D field
#-------------------------------------------------------------------------------

def get_plane(FIELD,plane,coordinate):

    if   plane == "xz":
        FIELD = FIELD[coordinate,:,:]             # (x,z)
    elif plane == "xy":
        FIELD = FIELD[:,coordinate,:]             # (x,y)
    elif plane == "yz":
        FIELD = FIELD[:,:,coordinate].transpose() # (y,z) (transpose so z (not y) is plotted down)

    return FIELD

#-------------------------------------------------------------------------------
# Calculate the lattice viscosity, thermal diffusivity, gravity and
# physical time step. Update Ra_p and Pr if unspecified, or
# calculate Pr from the physical parameters if Pr was unspecified
# and ka_p was specified. And set Ra_p to Ra if Ra_p unspecified
#-------------------------------------------------------------------------------

def get_lattice_viscosity(state,pturb,Ra,Pr,g_p,beta_p,dT_p,L_p,eta_p,rho_p,ka_p,g_l,L_l,nu_l,optimize_nu,nu_l_min,dt_p_input,Ma_sim,u_p_geo,gamma_T,rheology,geometry,print_lattice_info=True,rank=0):

    if dt_p_input != "unspecified" and dt_p_input != "calculate":
        error = False
        dt_p_int, timescale_dt_p , error = timevalue("dt_p",dt_p_input,error,rank)

        if timescale_dt_p != "timesteps":
            error = False
        if timescale_dt_p != "timesteps":
            dt_p       = timevalue_to_secs(dt_p_int,timescale_dt_p)
            dt_p_secs  = "known"
    else:
        dt_p_secs  = "unknown"

    Ra_l = Ra
    Pr_l = Pr

    if state == "solid":
        Ra_p = Ra_l
        Pr_p = "unknown"
    if state == "magma":
        Ra_p = "unknown"
        Pr_p = Pr_l

    C_p  = beta_p*dT_p
    C_l  = C_p
    dx_p = L_p/L_l      # Lattice resolution in physical units

    if print_lattice_info and rank == 0:
        print("==================================================================",flush=True)
        print("Getting lattice parameters from the run specifications"            ,flush=True)
        print("==================================================================",flush=True)
        print("state   =   %s            was specified" %(state)                  ,flush=True)
        print("------------------------------------------------------------------",flush=True)
        print("Ra_l    = %8.4g     was specified" %(Ra_l)                         ,flush=True)
    if state == "solid":
        Ra_p = Ra_l
        if print_lattice_info and rank == 0:
            print("Ra_p    = %8.4g     was set to Ra_l" %(Ra_p)                   ,flush=True)
            print("Pr_l    = %8.4g     was specified"   %(Pr_l)                   ,flush=True)
            print("------------------------------------------------------------------",flush=True)
            print("dt_p    = diffusion timescale = dt^d_p = L_l**2*ka_l/(L_l**2*ka_p)",flush=True)
            print("          provides the relationship to calculate ka_l",flush=True)
            print("------------------------------------------------------------------",flush=True)
    elif state == "magma":
        Pr_p = Pr_l
        if print_lattice_info and rank == 0:
            print("Pr_l    = %8.4g     was specified"   %(Pr_l),flush=True)
            print("Pr_p    = %8.4g     was set (Pr_p = Pr_l)" %(Pr_p),flush=True)
            print("------------------------------------------------------------------",flush=True)
        if isinstance(pturb,str) and pturb == "default":
            pturb_default = 3/7
            pturb         = pturb_default
            if print_lattice_info and rank == 0:
                print("pturb   = default      was specified",flush=True)
                print("p       = %8.4g     was defaulted to 3/7" %(pturb_default)     ,flush=True)
        else:
            print("p       = %8.4g     was specified" %(pturb),flush=True)
        if print_lattice_info and rank == 0:
            print("------------------------------------------------------------------",flush=True)
            print("dt^F_p  = Froude timescale = sqrt(g_l*C_l*L_p/(C_p*g_p*L_l)",flush=True)
            print("          provides the relationship to calculate g_l",flush=True)
            print("------------------------------------------------------------------",flush=True)
            print("dt_p    = the timescale of the simulation = dt^F_p*(Ra_p/Ra_l)**(1/2-p) due",flush=True)
            print("          to the Re propro Ra^p scaling at high Ra with p ~ 3/7 at Pr ~ 1",flush=True)
            print("------------------------------------------------------------------",flush=True)

    if print_lattice_info and rank == 0:
        len_geometry = len(geometry)
        print("geometry = %s" %(geometry),end="",flush=True)
        for i in range(12-len_geometry): print(" ",end="",flush=True)
        print("was specified",flush=True)
        print("------------------------------------------------------------------"    ,flush=True)
        print("L_p      = %8.4g    was specified" %(L_p),flush=True)
        if geometry == "cartesian":
            print("L_l      = %8.4g     was calculated for an %s geometry (L_l=nz-1)" %(L_l,geometry),flush=True)
        else:
            print("L_l      = %8.4g     was calculated for a %s geometry" %(L_l,geometry),flush=True)
        print("------------------------------------------------------------------"    ,flush=True)
        print("dx_p     = %8.4g     was calculated (dx_p=L_p/L_l)" %(dx_p)            ,flush=True)
        print("------------------------------------------------------------------"    ,flush=True)

        if isinstance(dt_p_input,float):
            print("dt_p     = %8.4g    was specified" %(dt_p_input),flush=True)
        else:
            len_dt_p = len(dt_p_input)
            print("dt_p    = %s" %(dt_p_input),end="",flush=True)
            for i in range(13-len_dt_p): print(" ",end="",flush=True)
            print("was specified",flush=True)
        if isinstance(Ma_sim,float):
            print("Ma_sim  = %8.4g     was specified" %(Ma_sim),flush=True)
        else:
            len_Ma_sim = len(Ma_sim)
            print("Ma_sim  = %s" %(Ma_sim),end="",flush=True)
            for i in range(13-len_Ma_sim): print(" ",end="",flush=True)
            print("was specified",flush=True)
        if isinstance(u_p_geo,float):
            print("u_p_geo = %8.4g     was specified" %(u_p_geo),flush=True)
        else:
            len_u_p_geo = len(u_p_geo)
            print("u_p_geo = %s" %(u_p_geo),end="",flush=True)
            for i in range(13-len_u_p_geo): print(" ",end="",flush=True)
            print("was specified",flush=True)
        print("------------------------------------------------------------------",flush=True)
        if isinstance(nu_l,float):
            print("nu_l    = %8.4g     was specified" %(nu_l)                     ,flush=True)
        else:
            len_nu_l = len(nu_l)
            print("nu_l    = %s" %(nu_l)                                   ,end="",flush=True)
            for i in range(13-len_nu_l): print(" "                         ,end="",flush=True)
            print("was specified"                                                 ,flush=True)
        if isinstance(g_l,float):
                print("g_l     = %8.4g     was specified" %(g_l)                  ,flush=True)
        else:
            len_g_l = len(g_l)
            print("g_l     = %s" %(g_l),end="",flush=True)
            for i in range(13-len_g_l): print(" "                          ,end="",flush=True)
            print("was specified"                                                 ,flush=True)
        if isinstance(eta_p,float):
            print("eta_p   = %8.4g     was specified" %(eta_p)                    ,flush=True)
        else:
            len_eta_p = len(eta_p)
            print("eta_p   = %s" %(eta_p)                                  ,end="",flush=True)
            for i in range(13-len_eta_p): print(" "                        ,end="",flush=True)
            print("was specified",flush=True)
        print("rho_p   = %8.4g     was specified" %(rho_p)                        ,flush=True)
        print("------------------------------------------------------------------",flush=True)
        print("g_p     = %8.4g     was specified" %(g_p)                          ,flush=True)
        print("ka_p    = %8.4g     was specified" %(ka_p)                         ,flush=True)
        print("------------------------------------------------------------------",flush=True)
        print("beta_p  = %8.4g     was specified" %(beta_p)                       ,flush=True)
        print("dT_p    = %8.4g     was specified" %(dT_p)                         ,flush=True)
        print("C_p     = %8.4g     was calculated (C_p=beta_p*dT_p)" %(C_p)       ,flush=True)
        print("==================================================================",flush=True)
    nu_p  = "unknown"
    if isinstance(eta_p,float):
        nu_p = eta_p/rho_p
        Pr_p = nu_p/ka_p
        if print_lattice_info and rank == 0:
            print("nu_p    = %8.4g    was calculated (nu_p=eta_p/rho_p)" %(nu_p)      ,flush=True)
            print("Pr_p    = %8.4g    was calculated (Pr_p=nu_p/ka_p)" %(Pr_p)        ,flush=True)
            print("------------------------------------------------------------------"    ,flush=True)
    if state == "solid" and Ra_p != "unknown" and np.abs((Ra_p-Ra_l)/Ra_l)>1.e-4:
        factor = np.sqrt(Ra_l/Ra_p)
        nu_p  *= factor
        Ra_p   = Ra_l
        Ra_p_check = g_p*C_p*dT_p*L_p**3/(nu_p*ka_p)
        if print_lattice_info and rank == 0:
            print("State   = %s, WARNING: Ra_p = %8.4g != Ra_l = %8.4g" %(state,Ra_p,Ra_l),flush=True)
            print("nu_p    = %8.4g    was re-calculated (nu_p *= sqrt(Ra_l/Ra_p)" %(nu_p) ,flush=True)
            print("Ra_p    = %8.4g    was re-set to Ra_l = %8.4g =? %8.4g = g_p*C_p*dT_p*L_p**3/(nu_p*ka_p)" %(Ra_p,Ra_p_check),flush=True)
            print("------------------------------------------------------------------"    ,flush=True)
    if state == "magma":
        if nu_p == "unknown":
            nu_p = ka_p * Pr_p
            if print_lattice_info and rank == 0: print("nu_p    = %8.4g    was calculated (ka_p * Pr_p)" %(nu_p),flush=True)
        Ra_p   = g_p*C_p*dT_p*L_p**3/(nu_p*ka_p)
        if print_lattice_info and rank == 0:
            print("Ra_p    = %8.4g   was calculated (g_p*C_p*dT_p*L_p**3/(nu_p*ka_p))" %(Ra_p),flush=True)
            print("------------------------------------------------------------------",flush=True)

#   Calculate the lattice parameters and physical timestep

    if state == "solid":
        if dt_p_secs == "known":
            ka_l = dt_p*L_l**2*ka_p/L_p**2
            nu_l = ka_l*Pr_l
            g_l  = nu_l**2*Ra_l/(C_l*L_l**3*Pr_l)
            if print_lattice_info and rank == 0:
                print("g_l     = %8.4g = nu_l**2*Ra_l/(C_l*L_l**3*Pr_l) was calculated" %(g_l) ,flush=True)
                print("nu_l    = %8.4g = ka_l*Pr_l                      was calculated" %(nu_l),flush=True)
                print("ka_l    = %8.4g = dt_p*L_l**2*ka_p/L_p**2        was calculated" %(ka_l),flush=True)

        elif dt_p_input == "calculate" and ((Ma_sim == "default" or Ma_sim != "unspecified") and (u_p_geo == "default" or u_p_geo != "unspecified")):
            if Ma_sim == "default":
                Ma_sim = 0.01  # Default to 0.01 for high accuracy
                if print_lattice_info and rank == 0: print("Ma_sim  = %8.3g cm/y was set to the default" %(Ma_sim),flush=True)
            else:
                if print_lattice_info and rank == 0: print("Ma_sim = %8.3g cm/y was specified" %(Ma_sim),flush=True)
            if u_p_geo   == "default":
                u_p_geo = 10.0 # Default to expected maximum velocity for solid state convection in cm/y 
                if print_lattice_info and rank == 0: print("u_p_geo = %8.3g cm/y was set to the default" %(u_p_geo),flush=True)
            else:
                if print_lattice_info and rank == 0: print("u_p_geo = %8.3g cm/y was specified" %(u_p_geo),flush=True)
            if print_lattice_info and rank == 0:
                print("------------------------------------------------------------------",flush=True)

            dt_p = Ma_sim * dx_p/((u_p_geo/u_cm_py) * np.sqrt(3))
            ka_l = dt_p*L_l**2*ka_p/L_p**2
            nu_l = ka_l*Pr_l
            g_l  = nu_l**2*Ra_l/(C_l*L_l**3*Pr_l)
            if print_lattice_info and rank == 0:
                print("dt_p    = %8.4g = Ma_sim * dx_p/((u_p_geo/u_cm_py)*sqrt(3)) was calculated" %(dt_p),flush=True) 
                print("ka_l    = %8.4g = dt_p*L_l**2*ka_p/L_p**2                   was calculated" %(ka_l),flush=True)
                print("nu_l    = %8.4g = ka_l*Pr_l                                 was calculated" %(nu_l),flush=True)
                print("g_l     = %8.4g = nu_l**2*Ra_l/(C_l*L_l**3*Pr_l)            was calculated" %(g_l) ,flush=True)
        elif nu_l != "unspecified":
            ka_l = nu_l/Pr_l
            g_l  = nu_l**2*Ra_l/(C_l*L_l**3*Pr_l)
            dt_p = L_p**2*ka_l/(L_l**2*ka_p)
            if print_lattice_info and rank == 0:
                print("nu_l    = %8.4g                                     was specified" %(nu_l) ,flush=True) 
                print("ka_l    = %8.4g = nu_l/Pr_l                         was calculated" %(ka_l),flush=True)
                print("g_l     = %8.4g = nu_l**2*Ra_l/(C_l*L_l**3*Pr_l)    was calculated" %(g_l) ,flush=True)
                print("dt_p    = %8.4g = L_p**2*ka_l/(L_l**2*ka_p)         was calculated" %(dt_p),flush=True)
        elif g_l != "unspecified":
            nu_l = np.sqrt(g_l*C_l*L_l**3*Pr_l/Ra_l)
            ka_l = nu_l/Pr_l
            dt_p = L_p**2*ka_l/(L_l**2*ka_p)
            if print_lattice_info and rank == 0:
                print("g_l     = %8.4g                                     was specified" %(g_l)  ,flush=True) 
                print("nu_l    = %8.4g = sqrt(g_l*C_l*L_l**3*Pr_l/Ra_l)    was calculated" %(nu_l),flush=True)
                print("ka_l    = %8.4g = nu_l/Pr_l                         was calculated" %(ka_l),flush=True)
                print("dt_p    = %8.4g = L_p**2*ka_l/(L_l**2*ka_p)         was calculated" %(dt_p),flush=True)

        dt_d_p = dt_p # Set the diffusion timescale

    elif state == "magma":
        if dt_p_secs == "known":
            g_l  = dt_p**2*g_p*C_p*L_l/(C_l*L_p)
            nu_l = np.sqrt(g_l*C_l*L_l**3*Pr_l/Ra_l)
            ka_l = nu_l/Pr_l
            if print_lattice_info and rank == 0:
                print("g_l     = %8.4g = dt^F_p**2*g_p*C_p*L_l/(C_l*L_p)   was calculated" %(g_l) ,flush=True)
                print("nu_l    = %8.4g = sqrt(g_l*C_l*L_l**3*Pr_l/Ra_l)    was calculated" %(nu_l),flush=True)
                print("ka_l    = %8.4g = nu_l/Pr_l                         was calculated" %(ka_l),flush=True)
                print("------------------------------------------------------------------"        ,flush=True)
        elif dt_p_input == "calculate" and ((Ma_sim == "default" or Ma_sim != "unspecified") and (u_p_geo == "default" or u_p_geo != "unspecified")):
            if Ma_sim == "default":
                Ma_sim = 0.01  # Default to 0.01 for high accuracy
                if print_lattice_info and rank == 0: print("Ma_sim  = %8.3g cm/y was set to the default" %(Ma_sim),flush=True)
            else:
                if print_lattice_info and rank == 0: print("Ma_sim = %8.3g cm/y was specified" %(Ma_sim),flush=True)
            if u_p_geo   == "default":
                u_p_geo = 35.0 # Default to expected maximum velocity for magma ocean dynamics m/s 
                if print_lattice_info and rank == 0: print("u_p_geo = %8.3g m/s was set to the default" %(u_p_geo),flush=True)
            else:
                if print_lattice_info and rank == 0: print("u_p_geo = %8.3g m/s was specified" %(u_p_geo),flush=True)
            if print_lattice_info and rank == 0:
                print("------------------------------------------------------------------",flush=True)

            dt_p = Ma_sim * dx_p/(u_p_geo * np.sqrt(3))
            g_l  = dt_p**2*g_p*C_p*L_l/(C_l*L_p)
            nu_l = np.sqrt(g_l*C_l*L_l**3*Pr_l/Ra_l)
            ka_l = nu_l/Pr_l
            if print_lattice_info and rank == 0:
                print("dt^F_p  = %8.4g = Ma_sim * dx_p/(u_p_geo * sqrt(3)) was calculated" %(dt_p),flush=True) 
                print("g_l     = %8.4g = dt^F_p**2*g_p*C_p*L_l/(C_l*L_p)   was calculated" %(g_l ),flush=True)
                print("nu_l    = %8.4g = sqrt(g_l*C_l*L_l**3*Pr_l/Ra_l)    was calculated" %(nu_l),flush=True)
                print("ka_l    = %8.4g = ka_l = nu_l/Pr_l                  was calculated" %(ka_l),flush=True)
        elif nu_l != "unspecified":
            g_l  = nu_l**2*Ra_l/(C_l*L_l**3*Pr_l)
            dt_p = np.sqrt(g_l*C_l*L_p/(g_p*C_p*L_l))
            ka_l = nu_l/Pr_l
            if print_lattice_info and rank == 0:
                print("nu_l    = %8.4g                                     was specified"  %(nu_l),flush=True)
                print("g_l     = %8.4g = nu_l**2*Ra_l/(C_l*L_l**3*Pr_l     was calculated" %(g_l ),flush=True)
                print("ka_l    = %8.4g = nu_l/Pr_l                         was calculated" %(ka_l),flush=True)
                print("dt^F_p  = %8.4g = sqrt(g_l*C_l*L_p/(g_p*C_p*L_l))   was calculated" %(dt_p),flush=True)
        elif g_l != "unspecified":
            nu_l = np.sqrt(g_l*C_l*L_l**3*Pr_l/Ra_l)
            dt_p = np.sqrt(g_l*C_l*L_p/(g_p*C_p*L_l))
            ka_l = nu_l/Pr_l
            if print_lattice_info and rank == 0:
                print("g_l     = %8.4g                                     was specified" %(g_l ),flush=True)
                print("nu_l    = %8.4g = sqrt(g_l*C_l*L_l**3*Pr_l/Ra_l)    was calculated" %(nu_l),flush=True)
                print("ka_l    = %8.4g = nu_l/Pr_l                         was calculated" %(ka_l),flush=True)
                print("dt^F_p  = %8.4g = sqrt(g_l*C_l*L_p/(g_p*C_p*L_l)    was calculated" %(dt_p),flush=True)

        if print_lattice_info and rank == 0:
            print("------------------------------------------------------------------"            ,flush=True)
         
        dt_F_p = dt_p # Set Froude the timescale
        if type(pturb) == float:
            correction = np.power(Ra_p/Ra_l,0.5-pturb)
            dt_p      *= correction
            if print_lattice_info and rank == 0:
                print("dt_p    = %8.4g                                     was calculated (dt^F_p*(Ra_p/Ra_l)**(1/2-p)" %(dt_p),flush=True)
                print("------------------------------------------------------------------"        ,flush=True)

        if eta_p == "unspecified":
            eta_p = rho_p * nu_p
            if print_lattice_info and rank == 0:
                print("eta_p   = %8.4g = rho_p * nu_p                      was calculated"%(eta_p),flush=True)

        if rheology != "linear":
            if print_lattice_info and rank == 0:
                print("------------------------------------------------------------------"        ,flush=True)
                print("Rheology != linear, stress units = F/A = (m*a)/A = (rho.dx**3*a)/dx**2 = rho*dx*a",flush=True)
                print("=> sigma_l  = sigma_p/(rho_p*dx_p*g_p) will be used to convert from physical to lattice stress (eg. Tau_y_l)",flush=True)

    if print_lattice_info and rank == 0 and dt_p_input != "unspecified":
        if isinstance(Ma_sim,float):
            mach_l    = Ma_sim
        else:
            mach_l    = 0.01
        print("------------------------------------------------------------------"                ,flush=True)
        v_l       = mach_l/np.sqrt(3)
        v_l       = mach_l/np.sqrt(3)
        v_p       = v_l*dx_p/dt_p
        v_p_cm_py = v_p * u_cm_py
        if state == "solid":
            print("Ma = v_l/sqrt(3) = %5.3f => v_p = %8.3g cm/y = v_l*dx_p/dt_p*secs/y*100" %(mach_l,v_p_cm_py),flush=True)
        else:
            print("Ma = v_l/sqrt(3) = %5.3f => v_p = %8.3g m/s = v_l*dx_p/dt_p" %(mach_l,v_p)   ,flush=True)

#   Print a message on acuracy of LBM with SRT if it may be in question

    if (print_lattice_info) and rank == 0 and dt_p_input != "unspecified":
        if type(optimize_nu) == bool:
            optimize_nu_l = optimize_nu
            if optimize_nu: nu_l_ideal = 0.5/3
        elif type(optimize_nu) == float:
            optimize_nu_l = True
            nu_l_ideal = optimize_nu
        if optimize_nu_l:

            nu_l0  = nu_l
            ka_l0  = ka_l
            g_l0   = g_l
            tau_l0 = 3*nu_l0 + 0.5

            factor = nu_l_ideal/nu_l
            nu_l  *= factor
            ka_l  *= factor
            g_l   *= factor
            tau_l  = 3*nu_l + 0.5

            if print_lattice_info and rank == 0:
                print("====================================================",flush=True)
                print("    optimize_nu=True => updating nu_l for optimal accuracy",flush=True)
                print("====================================================",flush=True)
                print("    factor = %8.3g" %(factor),flush=True)
                print("----------------------------------------------------",flush=True)
                print("    nu_l   = %8.3g -> %8.3g" %(nu_l0,nu_l),flush=True)
                print("    ka_l   = %8.3g -> %8.3g" %(ka_l0,ka_l),flush=True)
                print("    g_l    = %8.3g -> %8.3g" %(g_l0,g_l),flush=True)
                print("----------------------------------------------------",flush=True)
                print("    tau_l  = %8.3g -> %8.3g" %(tau_l0,tau_l),flush=True)

#       The nu_l sweet spot for accuracy of LBM with SRT is nu_l E [0.0167,.167] == tau_l E [0.55,1]
#       Experience suggests nu_l >= 0.0001 is likely stable and nu_l < 0.00005 is likely unstable

        tau_l = 3*nu_l + 0.5
        tau_l_sweet_min  = 0.55 ; nu_l_sweet_min = (tau_l_sweet_min-0.5)/3
        tau_l_sweet_max  = 1.00 ; nu_l_sweet_max = (tau_l_sweet_max-0.5)/3
        nu_l_max         = 1.00
        nu_l_min1_stable = 1e-4
        nu_l_min2_stable = 5e-5

        print("=================================================================",flush=True)
        print("Guidelines on LBM accuracy and how to optimize run performance"    ,flush=True)
        print("------------------------------------------------------------------",flush=True)
        print("1. Ma = u_l/c_s = u_l/sqrt(3) should be less than about 0.05 - 0.1",flush=True)
        print("   so compressible errors which scale as O(Ma^2) remain small and" ,flush=True)
        print("   do not swamp viscous accuracy."                                 ,flush=True)
        print("2. This is not a sharp limit, so for some problems, Ma = 0.1 may"  ,flush=True)
        print("   be accurate, but other problems may require lower Ma such as"   ,flush=True)
        print("   0.01 or less. If unsure, tests can be done at different Ma to"  ,flush=True)
        print("   determine how high Ma may be before results are impacted."      ,flush=True)
        print("3. Running at the highest accurate Ma minimizes run time. For"     ,flush=True)
        print("   example, runs at Ma = 0.1 is 10 x faster than at Ma = 0.01."    ,flush=True)
        print("   Since Ma propto dt_p, Ma can be controlled by adjusting dt_p."  ,flush=True)
        print("4. The lattice viscosity nu_l also affects accuracy. The accuracy" ,flush=True)
        print("   sweet spot is nu_l in the range 0.02 -> 0.2 as this maintains"  ,flush=True)
        print("   2-nd order spatial accuracy & avoids excessive numerical diffusion",flush=True)
        print("   while keeping relaxation fast enough to properly track equilibrium.",flush=True)
        print("5. Lower lattice viscosities to 0.001 or even 0.0001 are Ok but"   ,flush=True)
        print("   eventually, if nu_l is too low, the LBM will become unstable."  ,flush=True)
        print("   Our experience is that instability occurs around 0.0001 -> 0.00005." ,flush=True)
        print("6. Higher nu_l such that tau_l = nu_l/c_s^2 + 0.5 >> 1 where c_s"  ,flush=True)
        print("   = 1/sqrt(3) increases numerical dissipation & degrades no-slip" ,flush=True)
        print("   boundary conditions. However, free-slip BC are not affected."   ,flush=True)
        print("7. Performance can strongly depend on chips and the cache blocking",flush=True)
        print("   size. It is advised to run tests at different nxi_min & nxi_max",flush=True)
        print("   to obtain the optimal cache size for a given problem & machine.",flush=True)
        print("   The default settings obtain the most square or cubic domain size",flush=True)
        print("   as this minimizes communication costs between cores in theory." ,flush=True)
        print("   However, the optimal domain size depends on problem size & number",flush=True)
        print("   of cores. If domains are too small, communication costs are high",flush=True)
        print("   & cache blocking inefficient so the run will be slow. Conversely,",flush=True)
        print("   if domains are too large, the bus can be a bottlenech for feeding",flush=True)
        print("   data to the cache & runs may become slow. Hence, it is important",flush=True)
        print("   to run on the right number of cores for a given problem size"   ,flush=True)
        print("   hits the domain size sweet spot for optimal performance per core.",flush=True)
        print("   This can make a large difference such as a factor of 4x !!!!"   ,flush=True)
        print("------------------------------------------------------------------",flush=True)
        print("nu_l = %8.3g" %(nu_l),end="",flush=True)
        if nu_l > nu_l_max:
            print(" >> sweet-spot range for LBM accuracy => accuracy will be degraded.",flush=True)
        elif nu_l > nu_l_sweet_max:
            print(" > sweet-spot range but not by much => accuracy is likely Ok." ,flush=True)
        elif nu_l >= nu_l_sweet_min:
            print(" is within the sweet-spot for LBM accuracy."                   ,flush=True)
        elif nu_l >= nu_l_min1_stable:
            print(" < sweet-spot range but > %8.3g so accuracy & stability likely Ok.",flush=True)
        elif nu_l >= nu_l_min2_stable:
            print("is below the sweet-spot range & below %8.3g so accuracy & stability at risk.",flush=True)
        else:
            print("is below %8.3g so stability is at high risk!"                  ,flush=True)

    return Ra_p, Pr_p, eta_p, nu_p, ka_p, Ra_l, Pr_l, g_l, nu_l, ka_l, dx_p, dt_p, pturb

#-------------------------------------------------------------------------------
# Convert a field to the desired units ensuring the output data has no overflows
#-------------------------------------------------------------------------------

def convert_units(FIELD,dtype_data,c0,c1):

    FIELD = np.float32(FIELD)

    FIELD = c0 + c1 * FIELD

    return FIELD

#----------------------------------------------------------------------------
# Get the factor to convert velocity to output units & set the colorbar label
#
# If output_units = lattice, conversion = 1 & cbar label = "lattice units"
# If output_units = physical, convert to cm/year if dt_p > secs/year
# and otherwise convert to m/sec, and set the cbar label appropriately
#----------------------------------------------------------------------------

def get_u_convert(Pr,Ra,dx_p,dt_p,output_units="physical",latex=False):

#   Set the conversion factor depending on the output units and dt_p

    if output_units == "lattice":

        u_convert = 1

        if not latex: cbar_label_u = "Lattice units"
        else:         cbar_label_u = r'Lattice units'

    elif dt_p > sec_y:  # cm/year scale if dt_p >  1 year

        u_convert = dx_p/dt_p*u_cm_py

        if not latex: cbar_label_u = "cm/year"
        else:         cbar_label_u = r'cm/year'

    else:               # m/s     scale if dt_p <= 1 year

        u_convert = dx_p/dt_p

        if not latex: cbar_label_u = "m/s"
        else:         cbar_label_u = r'm/s'

    return u_convert, cbar_label_u

#===============================================================================
# Routines to get R & L from symbolic names and to initialize an annulus lattice
#===============================================================================

#-------------------------------------------------------------------------------
# Function to set R and L from descriptive names
#-------------------------------------------------------------------------------

def get_R_and_L(R_input,L_input,fouts,print_R_and_L_info=True,rank=0):

    if print_R_and_L_info and rank == 0 and (type(R_input) == str or type(L_input) == str):
        print("Getting R_p and/or L_p from specified symbolic names ...",flush=True)
        print("==================================================================",flush=True)

#   Set the radius of the planet/exoplanet/moon

    if type(R_input) == str:
        if   R_input == "Mercury":
            R = 2.437e6
        elif R_input == "Venus":
            R = 6.052e6
        elif R_input == "Earth":
            R = 6.378e6    # This is the equatorial radius. (cf. volumetric mean radius = 6371 and polar radius = 6356)
        elif R_input == "Moon":
            R = 1.737e6
        elif R_input == "Mars":
            R = 3.389e6
        elif R_input == "Uranus":
            R = 2.5362e7
        elif R_input == "Neptune":
            R = 2.4622e7
        else:
            if print_R_and_L_info and rank == 0:
                for i in range(len(fouts)):
                    print("ERROR: R = %s is unknown" %(R_input),file=fouts[i],flush=True)
            exit()
        if print_R_and_L_info and rank == 0:
            print("R = %8.4g m was set for \"%s\"" %(R,R_input),flush=True)

#   Set the thickness of the layer

    if type(L_input) == str:
        if type(R_input) == str:
            if R_input == "Earth" and (L_input == "mantle" or L_input == "upper_mantle"):
                if   L_input == "mantle"      : L = 2.9e6
                elif L_input == "upper_mantle": L = 6.6e5
                if print_R_and_L_info and rank == 0:
                    print("L = %8.4g m was set for \"%s\"" %(L,L_input),flush=True)
            else:
                if rank == 0:
                    for i in range(len(fouts)):
                        print("ERROR: L = %s is invalid for R = %s" %(L_input,R_input),file=sys.stderr,flush=True)
                exit()

    return R,L

#-------------------------------------------------------------------------------
# Function to get the annular lattice
#-------------------------------------------------------------------------------

def get_annular_lattice(R_p,L_p,nx,theta1,theta2,BC_sides,rank,print_annulus_info=True):

#   Initialize thetas array in radians with nx thetas at a regular theta increment dtheta from
#   theta1 to theta2 allowing for periodic side boundaries where f(thetas[0]) = f(theta[n*nx])
#   which means the left and right of the lattice are separated by dtheta. Hence, for
#   periodic sides, the thetas array will range from theta1+dtheta/2  to theta2-dtheta/2

    theta1_rad  = np.pi*theta1/180                  # Start theta  (radians)
    theta2_rad  = np.pi*theta2/180                  # End   theta  (radians)
    Dtheta_rad  = theta2_rad - theta1_rad           # Range theta  (radians)

    if BC_sides == "periodic":
        dtheta_rad  = Dtheta_rad/nx                 # dtheta       (radians) for peridic side boundaries
        theta1_rad += dtheta_rad/2                  # Update theta1_rad      for periodic side boundaries
        theta2_rad -= dtheta_rad/2                  # Update theta2_rad      for periodic side boundaries
    else:
        dtheta_rad  = Dtheta_rad/(nx-1)             # dtheta       (radians)
    thetas = np.linspace(theta1_rad,theta2_rad,nx)  # Thetas array (radians)
    theta1_deg  = 180*theta1_rad/np.pi              # Start theta  (degrees)
    theta2_deg  = 180*theta2_rad/np.pi              # End   theta  (degrees)
    dtheta_deg  = 180*dtheta_rad/np.pi              # dtheta       (degrees) 

#   Calculate dx_p, the distance on the surface corresponding to dtheta assuming
#   dtheta is small as required for the annulus to be well resolved by the TLBM

    X_p         = R_p*(theta2_rad-theta1_rad)       # Surface   extent    (m)
    dx_p        = X_p/(nx-1)                        # Lattice spacing at surface in m
    dx_p_km     = dx_p/1000                         # Lattice spacing at surface in km

#   Get the normalized radii at the surface and base of the layer and initialize
#   the transformation array xi spanning xi_0 at the base to 1 at the surface

    n           = 1000                              # Number of  xi indices spanning normalized depth z
    xi_0        = np.log((R_p-L_p)/R_p)/np.pi+1     # Initialize xi at the layer base
    xi          = np.linspace(xi_0,1,n)             # Initialize xi array
    r           = np.exp(np.pi*(xi-1))              # Normalized radii from the layer base to surface

#   Calculate the cubic spline of z(r) of the depth index as a function of normalized radius,
#   calculate the slope dz/dz of this spline, and calculate the lattice spacing dx(i) (ie. dx vs i)
#   where i is the normalized radius index ranging to n = number of points to resolve dx 

    z           = np.linspace(0,1,n)                # Normalized depth z where the base (z=0) is at i = 0 and surface (z=1) is at i = n-z
    cs_r_z      = CubicSpline(z,r)                  # Cubic spline of normalized radii r from base to surface
    dr_dz       = cs_r_z(z,1)                       # First derivative of r
    dx_rev      = dr_dz                             # Relative lattice spacing dx_rev[i] with dx_rev[n-1] at the surface
    dx          = np.zeros(n)
    for i in range(n):
        dx[i]   = dx_rev[n-1-i]                     # Reverse dx_rev[i] index order so dx[0] at surface and dx[n-1] at base
    dx          = dx/dx[0]                          # Scale dx so dx is unity at the surface and decreases with depth index i

#   Calculate the radius, layer thickness and surface extent in km, and
#   the surface, mid-depth and base cell sizes, extents and aspect ratios

    R_p_km      = R_p/1000                          # Radius              (km)
    L_p_km      = L_p/1000                          # Layer     thickness (km)
    dxu_l       = dx[0  ]                           # Surface   cell size (lu)
    dxd_l       = dx[n-1]                           # Base      cell size (lu)
    dxu_p_km    = dx_p_km*dx[0  ]                   # Surface   cell size (km)
    dxd_p_km    = dx_p_km*dx[n-1]                   # Base      cell size (km)
    Xu_p_km     = (nx-1)*dxu_p_km                   # Surface   extent    (km)
    Xd_p_km     = (nx-1)*dxd_p_km                   # Base      extent    (km)
    Xu_p_km     = (nx-1)*dxu_p_km                   # Surface   extent    (km)
    Xd_p_km     = (nx-1)*dxd_p_km                   # Base      extent    (km)
    Xm_p_km     = (Xu_p_km+Xd_p_km)/2               # Mid-depth extent    (km)
    Au          = Xu_p_km/L_p_km                    # Surface   aspect surface
    Ad          = Xd_p_km/L_p_km                    # Base      aspect base
    Am          = Xm_p_km/L_p_km                    # Mid-depth aspect ratio
    A           = Am                                # Aspect ratio = mid-depth aspect ratio

#   Calculate the cubic spline of dx and its slope dx/dz, integrate dx/dz to obtain L_l = L_p/dx_p (the lattice
#   layer thickness), and calculate the number of depths in the annular lattice nz (closest integer to L_l + 1)

    cs_dx       = CubicSpline(z,dx)                 # Cubic spline of dx
    int_dx      = cs_dx.integrate(0,1)              # Layer thickness in lattice units = integral (dx/dz) dz
    L_l         = (nx-1)*L_p_km/Xu_p_km             # Calculate L_l from the "upper" aspect ratio, so dx_p = L_p/L_l = dx_p at the surface
    nz_flt      = L_l/int_dx+1                      # Calculate nz as a float
    nz          = int(nz_flt+0.5)                   # Calculate nz for the annuluar lattice

#   Repeat the process to calculate dx(z) over nz, the now known lattice size

    xi          = np.linspace(xi_0,1,nz)            # Initialize xi array
    r           = np.exp(np.pi*(xi-1))              # Normalized radii from the layer base to surface
    z           = np.linspace(0,1,nz)               # Normalized depth z where the base is at iz = 0 and surface is at iz = nz-1
    cs_r_z      = CubicSpline(z,r)                  # Cubic spline of normalized radii r from base to surface
    dr_dz       = cs_r_z(z,1)                       # First derivative of r (z)
    dx_rev      = dr_dz                             # Lattice spacing dx_rev[z] with dx_rev[nz-1] at the surface
    dx          = np.zeros(nz)
    radii       = np.zeros(nz)
    for iz in range(nz):
        dx   [iz] = dx_rev[nz-1-iz]                 # Lattice spacing  with dx   [0] at the surface
        radii[iz] = r     [nz-1-iz]                 # Normalized radii with radii[0] at the surface
    dx          = dx/dx[0]                          # Scale dx so it is unity at surface
    radii      *= R_p                               # Radii in metres

#   Calculate the squareness measure S: lim(dtheta->0) S = 1 (ie. S ~ 1 for well resolved models)

    dx_diff     =  dx[0:nz-1] - dx[1:nz]            # Difference between adjacent cell edge sizes (dr       = cell size in radial direction midway between cells)
    dx_centre   = (dx[0:nz-1] + dx[1:nz])/2         # Cell centre sizes                           (r.dtheta = cell size in theta  direction midway between cells)
    dx_error    = dx_diff/dx_centre                 # Relative error between cell edge and cell centre sizes
    S           = np.max(1 + dx_error)              # Measure of closeness
    error_min   = np.min(dx_error)                  # Min error (ie. at the surface)
    error_max   = np.max(dx_error)                  # Max error (ie. at the base   )

    if print_annulus_info and rank == 0:

#       Plot the relative cell size vereus depth

        dx_l    = dx
        iz      = z*nz
        dx_p_km = dx_p*dx/1000
        z_p_km  = (R_p-radii)/1000

#       This is set to False as stdout sometimes becomes corrupted if it is True. ?Python issue? ?WS issue? ?Linux issue?

        plot_annulus_cell_size = False

        if plot_annulus_cell_size:
            plot_dir = "./"

            filename = plot_dir + "annulus_cell_size_vs_z-index_lattice.png"
            plt.title("Relative lattice cell size dx_p(z)/dx_p(0) vs depth index iz")
            plt.plot(iz,dx_l,"-",color="blue")
            plt.xlabel("Depth index")
            plt.ylabel("Relative lattice cell size")
            plt.savefig(filename)  ; plt.close()
            print("    Saved %s" %(filename),flush=True)

            filename = plot_dir + "annulus_cell_size_vs_z-index_physical.png"
            plt.title("Physical cell size dx_p vs depth index iz")
            plt.plot(iz,dx_p_km,"-",color="blue")
            plt.xlabel("Depth index")
            plt.ylabel("Cell size (km)")
            plt.savefig(filename)  ; plt.close()
            print("    Saved %s" %(filename),flush=True)

            filename = plot_dir + "annulus_cell_size_vs_z_lattice.png"
            plt.title("Relative lattice cell size dx_p(z)/dx_p(0) vs depth")
            plt.plot(z_p_km,dx_l,"-",color="blue")
            plt.xlabel("Depth (km)")
            plt.ylabel("Relative lattice cell size")
            plt.savefig(filename)  ; plt.close()
            print("    Saved %s" %(filename),flush=True)

            filename = plot_dir + "annulus_cell_size_vs_z_physical.png"
            plt.title("Physical cell size dx_p(z) vs depth")
            plt.plot(z_p_km,dx_p_km,"-",color="blue")
            plt.xlabel("Depth (km)")
            plt.ylabel("Cell size (km)")
            plt.savefig(filename)  ; plt.close()
            print("    Saved %s" %(filename),flush=True)

        print("=================================================================="    ,flush=True)
        print("R_p       = %7.2f km = radius of planet/expolanet"         %(R_p_km)   ,flush=True)
        print("------------------------------------------------------------------"    ,flush=True)
        print("L_p       = %7.2f km = annulus thickness"                  %(L_p_km)   ,flush=True)
        print("L_l       = %7.2f lu = annulus thickness in lattice units" %(L_l)      ,flush=True)
        print("------------------------------------------------------------------"    ,flush=True)
        print("nx        = %7d = number of thetas     in annulus"   %(nx)             ,flush=True)
        print("nz        = %7d = number of depths     in annulus"   %(nz)             ,flush=True)
        print("------------------------------------------------------------------"    ,flush=True)
        print("theta1    = %7.2f degrees = annulus start"           %(theta1)         ,flush=True)
        print("theta2    = %7.2f degrees = annulus end"             %(theta2)         ,flush=True)
        print("------------------------------------------------------------------"    ,flush=True)
        print("dxu_p     = %7.2f (km) = physical cell size at surface" %(dxu_p_km)    ,flush=True) 
        print("dxd_p     = %7.2f (km) = physical cell size at base"    %(dxd_p_km)    ,flush=True) 
        print("------------------------------------------------------------------"    ,flush=True)
        print("dxu_l     = %7.2f (lu) = lattice cell size at surface"       %(dxu_l)  ,flush=True) 
        print("dxd_l     = %7.2f (lu) = lattice cell size at base"          %(dxd_l)  ,flush=True) 
        print("------------------------------------------------------------------"    ,flush=True)
        print("Xu_p      = %7.1f km = annulus    extent at surface"       %(Xu_p_km)  ,flush=True)
        print("Xm_p      = %7.1f km = annulus    extent at mid-depth"     %(Xm_p_km)  ,flush=True)
        print("Xd_p      = %7.1f km = annulus    extent at base"          %(Xd_p_km)  ,flush=True)
        print("------------------------------------------------------------------"    ,flush=True)
        print("Au        = %7.2f = surface   aspect ratio = Xu_p/L_p"            %(Au),flush=True)
        print("Am        = %7.2f = mid-depth aspect ratio = Xm_p/L_p"            %(Am),flush=True)
        print("Ad        = %7.2f = base      aspect ratio = Xd_p/L_p"            %(Ad),flush=True)
        print("------------------------------------------------------------------"    ,flush=True)
        print("A         = %7.2f = aspect ratio = mid-depth aspect ratio"        %(A) ,flush=True)
        print("------------------------------------------------------------------"    ,flush=True)
        print("S         = %7.5f = cell squareness (dtheta small => well resolved model => S ~ 1)" %(S),flush=True)
        print("S_error   = %6.3f" %(100*(S - 1))                               ,end="",flush=True)
        print(" %"                                                                    ,flush=True)

    return np.float32(radii), np.float32(thetas), nz, np.float32(L_l),  np.float32(dx), np.float32(A)

