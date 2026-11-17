
import numpy as np
import sys

# Initialize conversion factors

sec_ns   = 1e-9            # Secs to nanosecs  conversion factor
sec_us   = 1e-6            # Secs to microsecs conversion factor
sec_ms   = 1e-3            # Secs to millisecs conversion factor
sec_m    = 60              # Seconds to Mins   conversion factor
sec_h    = 60    *sec_m    # Seconds to Hours  conversion factor
sec_d    = 24    *sec_h    # Seconds to Days   conversion factor
sec_w    = 7     *sec_d    # Seconds to Weeks  conversion factor
sec_y    = 365.25*sec_d    # Seconds to Years  conversion factor
sec_M    = sec_y/12        # Seconds to Months conversion factor
sec_ky   = 10**3 *sec_y    # Seconds to Ky     conversion factor
sec_my   = 10**6 *sec_y    # Seconds to My     conversion factor
sec_gy   = 10**9 *sec_y    # Seconds to Gy     conversion factor
u_cm_py  = 100   *sec_y    # m/sec to cm/y     conversion factor

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

#----------------------------------------------------------------------------
#   Get the time text to print from a time specification
#----------------------------------------------------------------------------

def get_time_value_text_int_timescale_from_time_spec(time_spec,time_value,dt_p):

    error = False
    time_int,time_scale,error    = timevalue(time_spec,time_value,error,0)
    time_steps,time_steps_float,time_steps_factor = timevalue_to_timesteps(time_int,time_scale,dt_p,0)

    return time_steps, time_int, time_scale

#-------------------------------------------------------------------------------
# Break down a timescale specification text into the integer & text parts
#-------------------------------------------------------------------------------

valid_timescales = ["timesteps", "Nsecs", "Usecs", "Msecs", "Secs", "Mins" ,"Hours", "Days", "Weeks", "Months", "Years", "Ky"  , "My"  , "Gy"  ]
factor_timesteps = [ 1         , sec_ns , sec_us , sec_ms , 1     , sec_m  , sec_h  , sec_d , sec_w , sec_M   , sec_y  , sec_ky, sec_my, sec_gy]
physical_timescales = valid_timescales                           [1:]
conversion_timesteps_to_physical_time = factor_timesteps         [1:]

#-------------------------------------------------------------------------------
# Convert an input time value in time scale units to timesteps
#-------------------------------------------------------------------------------

def timevalue_to_timesteps(value,timescale,dt_p,rank=0):

    if value != "unspecified":
        for i in range(len(valid_timescales)):
            if   timescale == "timesteps":
                timesteps = value
            elif timescale == valid_timescales[i]:

#               Add 1.5 to ensure last output is done, ie. at time step nt

                timesteps = int(value*factor_timesteps[i]/dt_p+1.5)
    else:
        timesteps = 0

    return timesteps

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
    error = True

    if len(parts) == 2 or len(parts) == 1:
        if len(parts) == 1:
            if timetext == "unspecified":
                timevalue = "unspecified"
                timescale = "unspecified"
                error     = False
            else:
                timevalue = int(parts[0]) ; timescale = "timesteps"
                error = False
        else:
            timevalue = "undefined"
            timescale = parts[1]
            for i in range(len(valid_timescales)):
                if timescale == valid_timescales[i]:
                    error = False
                    timevalue = int(parts[0])
        if error:
            for iout in range(2):
                if   iout == 0: output_print = sys.stderr
                elif iout == 1: output_print = sys.stdout
                if rank == 0:
                    for i_out in range(2):
                        if   i_out == 0: output_print = sys.stderr
                        elif i_out == 1: output_print = sys.stdout
                        print("ERROR: %s = %s is invalid. %s must be of form %s = n_timescale" %(parameter,timetext,parameter,timetext),file=output_print,flush=True)
                        print("       timescale is the time scale units (ie. timesteps, Secs, Hours, Days, Years,",file=output_print,flush=True)
                        print("       Ky, My or Gy) and n is an integer of how many units of timescale",file=output_print,flush=True)

    error = error_timescale or error

    return timevalue, timescale, error

#----------------------------------------------------------------------------
# Get seconds from time text in format n_timescale
#----------------------------------------------------------------------------

def timevalue_to_secs(timevalue,timescale,rank=0):

    if timescale != "unspecified":
        for i in range(len(valid_timescales)):
            if timescale == valid_timescales[i]:
                secs = factor_timesteps[i]*timevalue
    else:
        secs = timevalue

    return secs

#-------------------------------------------------------------------------------
# Get the text of the physical time for prints and plots
#
# Arguments are t (LBM time) and dt_p (physical time step)
#-------------------------------------------------------------------------------

def timesteps_to_physical_time_text(t,dt_p,output_units="physical"):

    if output_units == "lattice":
        time_text = "lattice units"

    t_p = t*dt_p          # Get the physical time in secs

    physical_time_text = time_to_physical_time_text(t_p,dt_p)

    return physical_time_text

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

def get_lattice_viscosity(state, power_correction,Ra,Pr,g_p,beta_p,dT_p,L_p,eta_p,rho_p,ka_p,Tau_y0_p,D,g_l,L_l,nu_l,nu_l_min,dt_p,gamma_T,rheology,rank):

    Ra_l = Ra
    Pr_l = Pr
    C_p  = beta_p*dT_p
    C_l  = C_p
    dx_p = L_p/L_l      # Lattice resolution in physical units
    if rank == 0:
        print("========================================================",flush=True)
        print("Calculating lattice parameters and/or dt_p as required",flush=True)
        print("========================================================",flush=True)
        if Ra_l != "unspecified":
            print("Ra_l    = %8.4g specified" %(Ra_l),flush=True)
        print("Pr_l    = %8.4g specified" %(Pr_l),flush=True)
        if eta_p != "unspecified":
            print("eta_p   = %8.4g specified" %(eta_p),flush=True)
        else:
            print("eta_p   = unspecified",flush=True)
        print("rho_p   = %8.4g specified" %(rho_p),flush=True)
        print("beta_p  = %8.4g specified" %(beta_p),flush=True)
        print("dT_p    = %8.4g specified" %(dT_p),flush=True)
        print("Setting C_l = C_p = beta_p*dT_p = %8.4g" %(C_p),flush=True)
        print("--------------------------------------------------------",flush=True)

    if eta_p != "unspecified":
        nu_p = eta_p/rho_p
        Ra_p = g_p*C_p*L_p**3/(nu_p*ka_p)
        Pr_p = nu_p/ka_p
        if rank == 0:
            print("All physical properties specified: eta_p, rho_p and ka_p",flush=True)
            print("    Calculated nu_p = %8.4g from eta_p and rho_p" %(nu_p),flush=True)
            print("    Calculated Ra_p = %8.4g from physical specifications" %(Ra_p),flush=True)
            print("    Calculated Pr_p = %8.4g from nu_p and ka_p" %(Pr_p),flush=True)

    if state == "solid":
        if isinstance(Ra,str) and Ra == "unspecified":
            if eta_p != "unspecified": Ra_l = Ra_p
        elif eta_p != "unspecified" and abs((Ra_p-Ra)/Ra)>0.001:
            eta_p0 = eta_p
            nu_p0  = nu_p
            factor = Ra_p/Ra
            eta_p *= factor
            nu_p  *= factor
            if rank == 0:
                print("State = solid WARNING: Ra_p = %8.4g != specified Ra = Ra_l = %8.4g" %(Ra_p,Ra_l),flush=True)
                print("    Correcting viscosity by (Ra_p/Ra_l) = %8.4g before setting Ra_p = Ra_l" %(factor),flush=True)
                print("    Specified eta_p = %8.4g corrected to eta_p  = %8.4g" %(eta_p0,eta_p),flush=True)
                print("    Specified nu_p  = %8.4g corrected to nu_p   = %8.4g" %(nu_p0 ,nu_p) ,flush=True)
        else:
            if rank == 0:
                if eta_p == "unspecified":
                    print("State = solid, eta_p = unspecified => setting Ra_p = Ra_l = %8.4g" %(Ra_l),flush=True)
                elif eta_p != "unspecified":
                    print("State = solid, eta_p = %8.3g => Ra_p = %8.3g = Ra_l" %(eta_p,Ra_p),flush=True)
        Ra_p = Ra_l
        Pr_p = g_p*C_p*L_p**3/(Ra_p*ka_p**2)
        nu_p = ka_p*Pr_p
        if rank == 0:
            print("    Calculated Pr_p = %8.4g from ka_p and Ra_p" %(Pr_p),flush=True)
            print("    Calculated nu_p  = %8.4g from ka_p  and Pr_p" %(nu_p),flush=True)

    elif state == "magma":

        if eta_p != "unspecified":
            Pr_p = nu_p/ka_p
            if abs((Pr_p-Pr_l)/Pr_l)>0.001:
                factor = Pr_l/Pr_p
                Pr_p  = Pr_l
                eta_p0 = eta_p
                nu_p0  = nu_p
                eta_p *= factor
                nu_p  *= factor
                if rank == 0:
                    print("State = magma, WARNING: Pr_p = %8.4g != specified Pr = Pr_l = %8.4g" %(Pr_p,Pr_l),flush=True)
                    print("    Correcting viscosity by to be consistent with Pr",flush=True)
                    print("    Specified eta_p = %8.4g corrected to eta_p  = %8.4g" %(eta_p0,eta_p),flush=True)
                    print("    Specified nu_p  = %8.4g corrected to nu_p   = %8.4g" %(nu_p0 ,nu_p) ,flush=True)
                    print("    Setting Pr_p = Pr_l = %8.4g" %(Pr_l),flush=True)
            else:
                if rank == 0:
                    print("State = magma, Pr_p = %8.4 = Pr_l = %8.4g" %(Pr_p,Pr_l),flush=True)
        else:
            Pr_p = Pr_l
            if rank == 0:
                print("State = magma => setting Pr_p = Pr_l = %8.4g" %(Pr_p),flush=True)

        Ra_p  = g_p*C_p*L_p**3/(Pr_p*ka_p**2)
        nu_p  = ka_p*Pr_p
        if rank == 0:
            print("    Calculated Ra_p  = %8.4g from ka_p  and Pr_p" %(Ra_p),flush=True)
            print("    Calculated nu_p  = %8.4g from ka_p  and Pr_p" %(nu_p),flush=True)

    eta_p = rho_p*nu_p
    if rank == 0:
        print("    Calculated eta_p = %8.4g from rho_p and nu_p" %(eta_p),flush=True)

#   Calculate the lattice parameters and physical timestep

    if state == "solid":
        if rank == 0:
            print("State = solid => using diffusion timescale: dt_p = L_p**2*ka_l/(L_l**2*ka_p)",flush=True)
        if dt_p != "unspecified":
            error = False ; rank = 0
            dt_p_int, timescale_dt_p , error = timevalue("dt_p",dt_p,error,rank)
            dt_p       = timevalue_to_secs(dt_p_int,timescale_dt_p)

            ka_l = dt_p*L_l**2*ka_p/L_p**2
            nu_l = ka_l*Pr_l
            g_l  = nu_l**2*Ra_l/(C_l*L_l**3*Pr_l)
            if rank == 0:
                print("    dt_p = %8.4g specified" %(dt_p),flush=True) 
                print("    Calculated ka_l = %8.4g" %(ka_l),flush=True)
                print("    Calculated nu_l = %8.4g" %(nu_l),flush=True)
                print("    Calculated g_l  = %8.4g" %(g_l),flush=True)
        elif nu_l != "unspecified":
            ka_l = nu_l/Pr_l
            g_l  = nu_l**2*Ra_l/(C_l*L_l**3*Pr_l)
            dt_p = L_p**2*ka_l/(L_l**2*ka_p)
            if rank == 0:
                print("    nu_l = %8.4g specified" %(nu_l),flush=True) 
                print("    Calculated ka_l = %8.4g" %(ka_l),flush=True)
                print("    Calculated g_l  = %8.4g" %(g_l),flush=True)
                print("    Calculated dt_p = %8.4g" %(dt_p),flush=True)
        elif g_l != "unspecified":
            nu_l = np.sqrt(g_l*C_l*L_l**3*Pr_l/Ra_l)
            ka_l = nu_l/Pr_l
            dt_p = L_p**2*ka_l/(L_l**2*ka_p)
            if rank == 0:
                print("    g_l = %8.4g specified" %(g_l),flush=True) 
                print("    Calculated nu_l = %8.4g" %(nu_l),flush=True)
                print("    Calculated ka_l = %8.4g" %(ka_l),flush=True)
                print("    Calculated dt_p = %8.4g" %(dt_p),flush=True)
    elif state == "magma":
        if rank == 0:
            print("State = magma => using Froude timescale: dt_p = sqrt(g_l*L_p/(g_p*L_l)",flush=True)
        if dt_p != "unspecified":
            g_l  = dt_p**2*g_p*C_p*L_l/(C_l*L_p)
            nu_l = np.sqrt(g_l*C_l*L_l**3*Pr_l/Ra_l)
            if rank == 0:
                print("dt_p = %8.4g specified" %(dt_p),flush= True)
                print("    Calculated g_l  = %8.4g" %(g_l),flush= True)
                print("    Calculated nu_l = %8.4g" %(nu_l),flush= True)
        elif nu_l != "unspecified":
            g_l  = nu_l**2*Ra_l/(C_l*L_l**3*Pr_l)
            dt_p = np.sqrt(g_l*C_l*L_p/(g_p*C_p*L_l))
            if rank == 0:
                print("nu_l = %8.4g specified" %(nu_l),flush= True)
                print("    Calculated g_l  = %8.4g" %(g_l),flush= True)
                print("    Calculated dt_p = %8.4g" %(dt_p),flush= True)
        elif g_l != "unspecified":
            nu_l = np.sqrt(g_l*C_l*L_l**3*Pr_l/Ra_l)
            dt_p = np.sqrt(g_l*C_l*L_p/(g_p*C_p*L_l))
            if rank == 0:
                print("g_l = %8.4g specified" %(g_l),flush= True)
                print("    Calculated nu_l = %8.4g" %(nu_l),flush= True)
                print("    Calculated dt_p = %8.4g" %(dt_p),flush= True)
        ka_l = nu_l/Pr_l
        if rank == 0:
            print("    Calculated ka_l = %8.4g" %(ka_l),flush= True)
        if power_correction != "unspecified" and abs(power_correction) > 1.e-6:
            dt_p0 = dt_p
            correction = np.power(Ra_p/Ra_l,power_correction)
            dt_p *= correction
            if rank == 0:
                print("--------------------------------------------------------",flush=True)
                print("power_correction   = %8.4g specified, dt_p0 = %8.4g will be corrected" %(power_correction,dt_p0),flush=True)
                print("    Corrected dt_p = %8.4g dt_p0*(Ra_p/Ra_l)**power_correction" %(dt_p),flush=True)
                print("    NOTE: This is an ad-hoc adjustment that can be used to tune magma ocean",flush=True)
                print("          velocities to be consistent with estimates. eg. 0.02 => ~ 50 m/s",flush=True)

#   Allow for nu_l_min oooOOO

    dt_dt0 = 1

    if rank == 0:
        print("========================================================",flush=True)
        print("Results: lattice parameters and dt_p",flush=True)
        print("========================================================",flush=True)
        print("dx_p = %8.3g m" %(dx_p),flush=True)
        print("dt_p = %8.4g s" %(dt_p),flush=True)
        print("nu_l = %8.4g, ka_l = %8.4g, g_l = %8.4g" %(nu_l,ka_l,g_l),flush=True)
        print("--------------------------------------------------------",flush=True)
        print("Ra_p = %8.4g, Ra_l = %8.4g" %(Ra_p,Ra_l),flush=True)
        print("Pr_p = %8.4g, Pr_l = %8.4g" %(Pr_p,Pr_l),flush=True)

#       Print out explanation of unit conversions

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
        print("dx_p                 = %8.3g m" %(dx_p),flush=True)
        print("dt_p                 = %8.3g s" %(dt_p),flush=True)
        print("------------------------------------------------------------------",flush=True)
        print("c0 = Pr_p/Pr_l       = %8.3g" %(c0),flush=True)
        print("c1 = sqrt(Ra_l/Ra_p) = %8.3g" %(c1),flush=True)
        if rheology != "linear":
            conversion_Tau_p_2_l  = (1/rho_p)
            conversion_Tau_p_2_l *= (g_l/g_p) / dx_p**(D-1)
            Tau_y0_l = conversion_Tau_p_2_l * Tau_y0_p
            print("------------------------------------------------------------------",flush=True)
            print("Stress units = F/A = m*a/dx**(D-1) where D = number of dimensions",flush=True)
            print("    stress    = (rho*dx**D)*dx/dt**2/dx**(D-1) = rho*dx**2/dt**2",flush=True)
            print("    rho_p                = %8.3g" %(rho_p),flush=True)
            print("    Tau_p_2_l = %8.3g = %8.3g = dt_p**2/(c0*rho_p*dx_p**2)" %(conversion_Tau_p_2_l,dt_p**2/(c0*rho_p*dx_p**2)),flush=True)
            print("    Tau_y0_p  = %8.3g" %(Tau_y0_p),flush=True)
            print("    Tau_y0_l  = %8.3g = Tau_p_2_l * Tau_y0_p" %(conversion_Tau_p_2_l * Tau_y0_p),flush=True)
        print("------------------------------------------------------------------",flush=True)
        print("g_p  = %8.3g = c0*g_l*dx_p/dt_p**2     where g_l  = %8.3g" %(c0*g_l*dx_p/dt_p**2,g_l),flush=True)
        print("------------------------------------------------------------------",flush=True)
        print("nu_p = %8.3g = c0*c1*nu_l*dx_p**2/dt_p where nu_l = %8.3g" %(nu_l*c0*c1*dx_p**2/dt_p,nu_l),flush=True)
        print("ka_p = %8.3g = c1*ka_l*dx_p**2/dt_p    where ka_l = %8.3g" %(c1*ka_l*dx_p**2/dt_p,ka_l),flush=True)
        print("------------------------------------------------------------------",flush=True)
        v_l       = 0.005
        v_p       = v_l*dx_p/dt_p
        v_p_cm_py = v_p * u_cm_py
        if state == "solid":
            print("v_p = v_l*dx_p/dt_p*secs/y*100 = %8.3g cm/y for v_l = %5.3f" %(v_p_cm_py,v_l),flush=True)
        else:
            print("v_p = v_l*dx_p/dt_p = %8.3g m/s for v_l = %5.3f" %(v_p,v_l),flush=True)

    return Ra_p, Pr_p, eta_p, nu_p, ka_p, Ra_l, Pr_l, g_l, nu_l, ka_l, dx_p, dt_p, dt_dt0

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

    if False: # oooOOO ooooooo

        if not latex:
            cbar_label_u = cbar_label_u + " (est. from Ra_p != Ra)"
        else:
            cbar_label_u = cbar_label_u + r' (est. from $Ra_p\ \neq\ Ra$)'

    return u_convert, cbar_label_u

#===============================================================================
# Routine to initialize an annular geometry
#===============================================================================

#-------------------------------------------------------------------------------
# Function to set R and L from descriptive names
#-------------------------------------------------------------------------------

def get_R_and_L(R,L):

#   Set the radius of the planet/exoplanet/moon

    if   R == "Mercury":
        R = 2.437e6
    elif R == "Venus":
        R = 6.052e6
    elif R == "Earth":
        R = 6.378e6
    elif R == "Moon":
        R = 1.737e6
    elif R == "Mars":
        R = 3.389e6
    elif R == "Uranus":
        R = 2.5362e7
    elif R == "Neptune":
        R = 2.4622e7

#   Set the thickness of the layer

    if   L == "whole_mantle": L = 2.9e6
    elif L == "upper_mantle": L = 6.6e5

#   Convert to float as necessary

    if isinstance(R,str): R = float(R)
    if isinstance(L,str): L = float(L)

    return R,L

#-------------------------------------------------------------------------------
# Function to get the annular geometry
#-------------------------------------------------------------------------------

def get_annular_geometry(R,L,theta1,theta2,n_theta,rank=0):

#   Get R and L values from names if necessary

    R, L = get_R_and_L(R,L)

#   Get the inner radius from the outer radius and layer thickness

    r_2 = R
    r_1 = r_2 - L

#   Convert from m to km

    convert = 1000
    r_1    /= convert
    r_2    /= convert
    R       = r_2
    L      /= convert

#   Set theta and d_theta

    theta1 = np.pi*theta1/180
    theta2 = np.pi*theta2/180
    d_theta = (theta2-theta1)/n_theta
    theta   = np.linspace(theta1,theta2-d_theta,n_theta)

#   Initialize xi and r for a layer of thickness one lattice spacing

    n_r  = 2                         # Initialize to layer thickness of one lattice spacing

    dr   = r_2*d_theta               # Initialize dr

    xi_1 = np.log(r_1/r_2)/np.pi+1   # Calculate  xi at base of layer (xi = 1 is at the top)
    xi   = np.linspace(xi_1,1,n_r)   # Initialize xi
    r    = r_2*np.exp(np.pi*(xi-1))  # Initialize r from xi
    dr   = r[1:n_r] - r[0:n_r-1]     # Initialize dr
    dx   = r[:n_r-1]*d_theta         # Initialize dx

    while  dr[n_r-2] > dx[n_r-2]:

        n_r += 1                         # Increment n_r
        squareness = np.zeros(n_r-1)

        xi   = np.linspace(xi_1,1,n_r)   # Reinitialize xi
        r    = r_2*np.exp(np.pi*(xi-1))  # Reinitialize r from xi
        dr   = r[1:n_r] - r[0:n_r-1]     # Recalculate dr
        dx   = r[:n_r-1]*d_theta         # Recalculate dx

#   Update the transformation using r_1 exactly

    n_r -= 1

    xi   = np.linspace(xi_1,1,n_r)   # Reinitialize xi
    r    = r_2*np.exp(np.pi*(xi-1))  # Reinitialize r from xi

#   Re-calculate dr as centered finite-difference where possible

    dr = np.zeros(n_r)
    for i_r in range(n_r):
        if i_r == 0:
            dr[i_r] =  r[i_r+1] - r[i_r]
        elif i_r == n_r-1:
            dr[i_r] =  r[i_r  ] - r[i_r-1]
        else:
            dr[i_r] = (r[i_r+1] - r[i_r-1])/2

#   Calculate the layer thickness in lattice units

    dx0 = r[n_r-1]*d_theta
    L_l = (r_2-r_1)/dx0

#   Calculate and print annular lattice information

    i_r = max(n_r//2,1)
    dx  = r[i_r]*d_theta
    dy  = dr[i_r]

    squareness = dx/dy

    ir_top  = n_r-1 ; r_top  = r[ir_top ] ; dx_top  = r[ir_top ]*d_theta
    ir_base = 0     ; r_base = r[ir_base] ; dx_base = r[ir_base]*d_theta

    if rank == 0:
        print("==================================================================",flush=True)
        print("Annular lattice summary",flush=True)
        print("==================================================================",flush=True)
        print("L_l     = %7.2f ,  n_r = %7d, n_theta = %d" %(L_l,n_r,n_theta),flush=True)
        print("------------------------------------------------------------------",flush=True)
        print("R       = %7d km,  L = %7d km" %(R,L),flush=True)
        print("------------------------------------------------------------------",flush=True)
        print("R_top   = %d km, dx_top  = %6.2g km" %(r_top ,dx_top ),flush=True)
        print("R_base  = %d km, dx_base = %6.2g km" %(r_base,dx_base),flush=True)
        print("------------------------------------------------------------------",flush=True)
        print("dx_base/dx_top = %5.3g, squareness = %6.4f" %(dx_base/dx_top,squareness),flush=True)
        print("------------------------------------------------------------------",flush=True)

#   Correct for surface being at i_r = n_r-1 in theabove rather than at i_r = 0

    dr     = np.float64(dr)
    dr_out = np.zeros(n_r,dtype=np.float64)
    r_out  = np.zeros(n_r,dtype=np.float64)
    for i_r in range(n_r):
        r_out [i_r] = r [n_r-i_r-1]
        dr_out[i_r] = dr[n_r-i_r-1]

#   Calculate factor to apply to lattice viscosity due to variable lattice spacing in r

    nu_factor = (dr_out[0]/dr_out)**2

#   Convert back to m from km

    r_out *= convert
    R     *= convert
    L     *= convert

    return L_l,n_r,R,L,r_out,theta,nu_factor,dr,r*d_theta
