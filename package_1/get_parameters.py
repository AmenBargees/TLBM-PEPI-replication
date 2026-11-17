
#-------------------------------------------------------------------------------
# Function to parameters from an input file with anything after a # being
# ignored and considered as a comment. The input file must have the format of
# parameter_name = parameter_value where parameter_value can be a string,
# float, int or bool represented by the string True of False.
#-------------------------------------------------------------------------------

import sys

# Function to replace underscores with spaces in a string

def str2str(value):

    output_value = []
    n = 0
    last_character = ""
    for i in range(len(value)):
        if value[i] == "_" and last_character != "\\":
            output_value.append(" ")
        elif value[i] != "\\":
            output_value.append(value[i])
        last_character = str(value[i])

        n += 1

    output_value = "".join(output_value)
    output_value = str(output_value)

    return output_value

# Function to convert a string True or False into a bool

def str2bool(value):

    global parameter, iline, line, nwarning

    if   value == "True"  or value == "true":
        result = True
    elif value == "False" or value == "false":
        result = False
    else:
        result = False
        nwarning += 1
        if rank == 0:
            for iout in range(2):
                if   iout == 0: output_print = sys.stderr
                elif iout == 1: output_print = sys.stdout
                print('WARNING: line number = %d, line = "%s"' %(iline,line.rstrip()),file=output_print,flush=True)
                print('         Must specify True or False for boolean variable "%s"' %(parameter),file=output_print,flush=True)
                print('         "%s" is being set to False' %(parameter),file=output_print,flush=True)
    return result

def get_parameters(filename,rank=0):

    global parameter, iline, line, nwarning

#   Run title and description for plots

    global run_title, description, description_color, description_box, description_location

#   Model size, time steps and geometry

    global L_p, nx, ny, nz, nt, geometry, theta1, theta2, n_TBL

#   Physical specifications

    global state, power_correction, Ra, Pr, Peclet, R_p, L_p, rho_p, eta_p, ka_p, g_p, beta_p, T0_p, dT_p
    global rho_z_file, P_z_file, T_z_file, L_z_file
    global dT_p_random, initialize_ncells, cell_amplitude, amp_sinusoid, x_periods, y_periods
    global H, Cp, temperature, gamma_T, rheology, nu_formula, Tau_y0_p, Tau_yi_p, nu_star, N
    global RKTLBM, collision, melt_time, melt_position, melt_extraction 
    global rho_p_TTG, rho1_p, rho_p_peridotite, rho2_p, rho_p_basalt, rho3_p, rho_r_file, rho_b_file, A, interface_thickness, nu_r_factor, beta_g_factor, g_comp_factor, w_diagonal, theta_u, theta_d, theta_l, theta_r
    global rayleigh_taylor_option, depth_light, perturbation_size, optimize_nu

#   Boundary condition specifications

    global BC_sides, BC_upper, BC_lower

#   Physical time step parameters
#       If g  is specified, nu will be calculated from g
#       If nu is specified, g  will be calculated from nu
#       If both nu and g are specified, nu will take precedence

    global time_per_timestep, nu_l, g_l, factor_dt

#   Mach control and accuracy/stability parameters

    global max_mach_limit, min_mach_limit, Inc_vary_dt, start_vary_dt
    global factor_increase_dt,factor_decrease_dt,factor_max_vary_dt,mach_limit_to_exit
    global nu_l_min

#   Accuracy parameters

    global float_bits, edge_bits, data_bits, shifted_distribution

#   Output specifiations

    global output_plots, output_data, output_timing, output_units

#   Output increment parameters

    global Inc_data, Inc_plots, Inc_bmk, Inc_timing, Inc_mach

#   Output plot parameters

    global plane, coordinates
    global plot_size, plot_size_factor, plot_colorbar, colormap, latex
    global plot_rho, plot_T, plot_u, plot_nu_total, plot_nu_kinematic, plot_nu_dynamic, plot_vectors, plot_Tau
    global plot_rho_eff, plot_L, plot_rho_s, plot_M_pd, plot_M_bs, plot_M_LP, plot_M_MP, plot_M_HP, plot_ce_LP, plot_ce_MP, plot_ce_HP, plot_ce_pd, plot_ce_bs
    global plot_T_sol_pd, plot_T_sol_bs, plot_T_liq_pd, plot_T_liq_bs, plot_P, plot_DI
    global plot_rho_r, plot_rho_b, plot_f_r, plot_f_b, plot_FF_mag, plot_FF_x, plot_FF_z, plot_FF_y
 
#   Output velocity vector and plot scale parameters

    global vector_overlay, vector_color, vector_length, vector_width, vector_interval, annulus_interval
    global clip_rho, gpow_rho, clip_T, gpow_T, fix_scale_T, clip_u, gpow_u
    global clip_nu_total, gpow_nu_total, clip_Tau, gpow_Tau

#   Performance parameters

    global code, overlap, domain_decomposition
    global nxi_min, nxi_max, nyi_min, nyi_max, nzi_min, nzi_max

#-------------------------------------------------------------------------------
#   Get read and set parameters
#-------------------------------------------------------------------------------

    if rank == 0:
        for iout in range(2):
            if   iout == 0: output_print = sys.stderr
            elif iout == 1: output_print = sys.stdout
            print("Reading input file ...",end="",file=output_print,flush=True)

    file_input  = open(filename,"r")
    
    iline    = 0       # Line number being read
    nerror   = 0       # Number of input errors
    nwarning = 0       # Number of input warnings
    verbose  = False   # Verbose mode prints out each parameter as it is set

    for line in file_input:
        iline += 1
        for i in range(len(line)):
            if line[i] == "#":        # Strip off comments
                line = line[0:i]
                break
        parameters = line.split()
        if (len(parameters) == 3 and parameters[1] != "=")  or \
           (len(parameters) != 3 and len(parameters) > 0):
                nerror += 1 
                if rank == 0:
                    for iout in range(2):
                        if   iout == 0: output_print = sys.stderr
                        elif iout == 1: output_print = sys.stdout
                        if nerror == 1: print()
                        print('ERROR: line number = %d, line = %s' %(iline,line.rstrip()),file=output_print,flush=True)
        elif len(parameters) == 3:
            parameter = parameters[0]
            value     = parameters[2]
            valid     = False
            if rank == 0:
                for iout in range(2):
                    if   iout == 0: output_print = sys.stderr
                    elif iout == 1: output_print = sys.stdout
                    if verbose:
                        if nerror == 1: print()
                        print("Setting parameter ",parameter," to ",value,file=output_print,flush=True)

#-------------------------------------------------------------------------------
#           Run title and description for plots
#-------------------------------------------------------------------------------

            if parameter == "run_title":
                valid = True ; run_title = value
            elif parameter == "description":
                valid = True ; description = value
            elif parameter == "description_color":
                valid = True ; description_color = str(value)
            elif parameter == "description_box":
                valid = True ; description_box = str2bool(value)
            elif parameter == "description_location":
                valid = True ; description_location = str(value)

#-------------------------------------------------------------------------------
#           Model size, time steps and geometry
#-------------------------------------------------------------------------------

            elif parameter == "R_p":
                valid = True ; R_p = str(value)
            elif parameter == "L_p":
                valid = True ; L_p = str(value)
            elif parameter == "nx":
                valid = True ; nx = int(value)
            elif parameter == "ny":
                valid = True ; ny = int(value)
            elif parameter == "nz":
                valid = True ; nz = int(value)
            elif parameter == "nt":
                valid = True ; nt = str(value)
            elif parameter == "geometry":
                valid = True ; geometry = str(value)
            elif parameter == "theta1":
                valid = True ; theta1 = float(value)
            elif parameter == "theta2":
                valid = True ; theta2 = float(value)
            elif parameter == "n_TBL":
                valid = True ; n_TBL = int(value)

#-------------------------------------------------------------------------------
#           Physical specifications
#-------------------------------------------------------------------------------

            elif parameter == "state":
                valid = True
                state = str(value)
            elif parameter == "power_correction":
                valid = True
                if value != "unspecified": power_correction = float(value)
                else:                      power_correction = str(value)
            elif parameter == "Ra":
                valid = True
                if value != "unspecified": Ra = float(value)
                else:                      Ra = str(value)
            elif parameter == "Pr":
                valid = True
                if value != "unspecified": Pr = float(value)
                else:                      Pr = str(value)
            elif parameter == "Peclet":
                valid = True ; Peclet = float(value)
            elif parameter == "rho_p":
                valid = True ; rho_p = float(value)
            elif parameter == "rho_z_file":
                valid = True ; rho_z_file = str(value)
            elif parameter == "P_z_file":
                valid = True ; P_z_file = str(value)
            elif parameter == "T_z_file":
                valid = True ; T_z_file = str(value)
            elif parameter == "L_z_file":
                valid = True ; L_z_file = str(value)
            elif parameter == "eta_p":
                valid = True
                if value != "unspecified":
                    eta_p = float(value)
                else:
                    eta_p = value
            elif parameter == "ka_p":
                valid = True
                if value != "unspecified": ka_p = float(value)
                else:                      ka_p = str(value)
            elif parameter == "g_p":
                valid = True ; g_p = float(value)
            elif parameter == "beta_p":
                valid = True ; beta_p = float(value)
            elif parameter == "T0_p":
                valid = True ; T0_p = float(value)
            elif parameter == "dT_p":
                valid = True ; dT_p = float(value)
            elif parameter == "dT_p_random":
                valid = True ; dT_p_random = float(value)
            elif parameter == "initialize_ncells":
                valid = True ; initialize_ncells = int(value)
            elif parameter == "cell_amplitude":
                valid = True ; cell_amplitude = float(value)
            elif parameter == "amp_sinusoid":
                valid = True ; amp_sinusoid = float(value)
            elif parameter == "x_periods":
                valid = True ; x_periods = float(value)
            elif parameter == "y_periods":
                valid = True ; y_periods = float(value)
            elif parameter == "H":
                valid = True ; H = float(value)
            elif parameter == "Cp":
                valid = True ; Cp = float(value)
            elif parameter == "temperature":
                valid = True ; temperature = str(value)
            elif parameter == "gamma_T":
                valid = True ; gamma_T = float(value)
            elif parameter == "rheology":
                valid = True ; rheology = str(value)
            elif parameter == "nu_formula":
                valid = True ; nu_formula = str(value)
            elif parameter == "Tau_y0_p":
                valid = True ; Tau_y0_p = float(value)
            elif parameter == "Tau_yi_p":
                valid = True ; Tau_yi_p = float(value)
            elif parameter == "nu_star":
                valid = True ; nu_star = float(value)
            elif parameter == "N":
                valid = True ; N = float(value)
            elif parameter == "RKTLBM":
                valid = True ; RKTLBM = str2bool(value)
            elif parameter == "collision":
                valid = True ; collision = str2bool(value)
            elif parameter == "melt_time":
                valid = True ; melt_time = float(value)
            elif parameter == "melt_position":
                valid = True ; melt_position = str(value)
            elif parameter == "melt_extraction":
                valid = True ; melt_extraction = str(value) 
            elif parameter == "rho_r_file":
                valid = True ; rho_r_file = str(value)
            elif parameter == "rho_b_file":
                valid = True ; rho_b_file = str(value)
            elif parameter == "rho_p_TTG" or parameter == "rho1_p":
                valid = True ; rho_p_TTG = float(value) ; rho1_p = rho_p_TTG
            elif parameter == "rho_p_peridotite" or parameter == "rho2_p":
                valid = True ; rho_p_peridotite = float(value) ; rho2_p = rho_p_peridotite
            elif parameter == "rho_p_basalt" or parameter == "rho3_p":
                valid = True ; rho_p_basalt = float(value) ; rho3_p = rho_p_basalt
            elif parameter == "A":
                valid = True ; A = float(value)
            elif parameter == "interface_thickness":
                valid = True ; interface_thickness = float(value)
            elif parameter == "nu_r_factor":
                valid = True ; nu_r_factor = float(value)
            elif parameter == "beta_g_factor": 
                valid = True ; beta_g_factor = float(value)
            elif parameter == "g_comp_factor": 
                valid = True ; g_comp_factor = float(value)
            elif parameter == "w_diagonal": 
                valid = True ; w_diagonal = float(value)
            elif parameter == "theta_u": 
                valid = True ; theta_u = float(value)
            elif parameter == "theta_d": 
                valid = True ; theta_d = float(value)
            elif parameter == "theta_l": 
                valid = True ; theta_l = float(value)
            elif parameter == "theta_r": 
                valid = True ; theta_r = float(value)
            elif parameter == "rayleigh_taylor_option": 
                valid = True ; rayleigh_taylor_option = str(value)
            elif parameter == "perturbation_size": 
                valid = True ; perturbation_size = float(value)
            elif parameter == "depth_light": 
                valid = True ; depth_light = float(value)
            elif parameter == "optimize_nu": 
                valid = True ; optimize_nu = str2bool(value)

#-------------------------------------------------------------------------------
#           Boundary condition specifications
#-------------------------------------------------------------------------------

            elif parameter == "BC_sides":
                valid = True ; BC_sides = str(value)
            elif parameter == "BC_upper":
                valid = True ; BC_upper = str(value)
            elif parameter == "BC_lower":
                valid = True ; BC_lower = str(value)

#-------------------------------------------------------------------------------
#           Physical time step specifications
#
#               If g  is specified, nu_f will be calculated from g
#               If nu is specified, g    will be calculated from nu
#               If both nu and g are specified, nu will take precedence
#               if time_per_timestep is specified, it takes overall precedence
#
#-------------------------------------------------------------------------------

            elif parameter == "time_per_timestep":
                valid = True
                time_per_timestep = str(value)
            elif parameter == "nu_l":
                valid = True
                if value != "unspecified": nu_l = float(value)
                else:                      nu_l = str(value)
            elif parameter == "g_l":
                valid = True
                if value != "unspecified": g_l = float(value)
                else:                      g_l = str(value)
            elif parameter == "factor_dt":
                valid = True ; factor_dt = float(value)

#-------------------------------------------------------------------------------
#           Mach number control specifications and minimum nu_l
#-------------------------------------------------------------------------------

            elif parameter == "start_vary_dt":
                valid = True ; start_vary_dt = str(value)
            elif parameter == "Inc_vary_dt":
                valid = True ; Inc_vary_dt = str(value)
            elif parameter == "max_mach_limit":
                valid = True ; max_mach_limit = float(value)
            elif parameter == "min_mach_limit":
                valid = True ; min_mach_limit = float(value)
            elif parameter == "factor_increase_dt":
                valid = True ; factor_increase_dt = float(value)
            elif parameter == "factor_decrease_dt":
                valid = True ; factor_decrease_dt = float(value)
            elif parameter == "factor_max_vary_dt":
                valid = True ; factor_max_vary_dt = float(value)
            elif parameter == "mach_limit_to_exit":
                valid = True ; mach_limit_to_exit = float(value)
            elif parameter == "nu_l_min":
                valid = True ; nu_l_min = float(value)

#-------------------------------------------------------------------------------
#           Accuracy specifications
#-------------------------------------------------------------------------------

            elif parameter == "float_bits":
                valid = True ; float_bits = int(value)
            elif parameter == "edge_bits":
                valid = True ; edge_bits = int(value)
            elif parameter == "data_bits":
                valid = True ; data_bits = int(value)
            elif parameter == "shifted_distribution":
                valid = True ; shifted_distribution = str2bool(value) ;

#-------------------------------------------------------------------------------
#           Output specifiations
#-------------------------------------------------------------------------------

            elif parameter == "output_plots":
                valid = True ; output_plots = str2bool(value)
            elif parameter == "output_data":
                valid = True ; output_data = str2bool(value)
            elif parameter == "output_timing":
                valid = True ; output_timing = str2bool(value)
            elif parameter == "output_units":
                valid = True ; output_units = str(value)

#-------------------------------------------------------------------------------
#           Output increment specifications
#-------------------------------------------------------------------------------

            elif parameter == "Inc_plots":
                valid = True ; Inc_plots = str(value)
            elif parameter == "Inc_bmk":
                valid = True ; Inc_bmk = str(value)
            elif parameter == "Inc_data":
                valid = True ; Inc_data = str(value)
            elif parameter == "Inc_mach":
                valid = True ; Inc_mach = str(value)
            elif parameter == "Inc_timing":
                valid = True ; Inc_timing = str(value)

#-------------------------------------------------------------------------------
#           Output plot specifications
#-------------------------------------------------------------------------------

            elif parameter == "plane":
                valid = True ; plane = str(value)
            elif parameter == "coordinates":
                valid = True ; coordinates = str(value)
            elif parameter == "plot_size":
                valid = True ; plot_size = str(value)
            elif parameter == "plot_size_factor":
                valid = True ; plot_size_factor = float(value)
            elif parameter == "plot_colorbar":
                valid = True ; plot_colorbar = str2bool(value)
            elif parameter == "colormap":
                valid = True ; colormap = str(value)
            elif parameter == "latex":
                valid = True ; latex = str2bool(value)
            elif parameter == "plot_rho":
                valid = True ; plot_rho = str2bool(value)
            elif parameter == "plot_T":
                valid = True ; plot_T = str2bool(value)
            elif parameter == "plot_u":
                valid = True ; plot_u = str2bool(value)
            elif parameter == "plot_nu_total":
                valid = True ; plot_nu_total = str2bool(value)
            elif parameter == "plot_nu_kinematic":
                valid = True ; plot_nu_kinematic = str2bool(value)
            elif parameter == "plot_nu_dynamic":
                valid = True ; plot_nu_dynamic = str2bool(value)
            elif parameter == "plot_vectors":
                valid = True ; plot_vectors = str2bool(value)
            elif parameter == "plot_Tau":
                valid = True ; plot_Tau = str2bool(value)
            elif parameter == "plot_L":
                valid = True ; plot_L = str2bool(value)
            elif parameter == "plot_rho_eff":
                valid = True ; plot_rho_eff = str2bool(value)
            elif parameter == "plot_rho_s":
                valid = True ; plot_rho_s = str2bool(value)
            elif parameter == "plot_M_pd":
                valid = True ; plot_M_pd = str2bool(value)
            elif parameter == "plot_M_bs":
                valid = True ; plot_M_bs = str2bool(value)
            elif parameter == "plot_M_LP":
                valid = True ; plot_M_LP = str2bool(value)
            elif parameter == "plot_M_MP":
                valid = True ; plot_M_MP = str2bool(value)
            elif parameter == "plot_M_HP":
                valid = True ; plot_M_HP = str2bool(value)
            elif parameter == "plot_ce_LP":
                 valid = True ; plot_ce_LP = str2bool(value)
            elif parameter == "plot_ce_MP":
                 valid = True ; plot_ce_MP = str2bool(value)
            elif parameter == "plot_ce_HP":
                 valid = True ; plot_ce_HP = str2bool(value)
            elif parameter == "plot_ce_pd":
                valid = True ; plot_ce_pd = str2bool(value)
            elif parameter == "plot_ce_bs":
                valid = True ; plot_ce_bs = str2bool(value)
            elif parameter == "plot_T_sol_pd":
                valid = True ; plot_T_sol_pd = str2bool(value)
            elif parameter == "plot_T_sol_bs":
                valid = True ; plot_T_sol_bs = str2bool(value)
            elif parameter == "plot_T_liq_pd":
                valid = True ; plot_T_liq_pd = str2bool(value)
            elif parameter == "plot_T_liq_bs":
                valid = True ; plot_T_liq_bs = str2bool(value)
            elif parameter == "plot_P":
                valid = True ; plot_P = str2bool(value)
            elif parameter == "plot_DI":
                valid = True ; plot_DI = str2bool(value)
            elif parameter == "plot_rho_r":
                valid = True ; plot_rho_r = str2bool(value)
            elif parameter == "plot_rho_b":
                valid = True ; plot_rho_b = str2bool(value)
            elif parameter == "plot_f_r":
                valid = True ; plot_f_r = str2bool(value)
            elif parameter == "plot_f_b":
                valid = True ; plot_f_b = str2bool(value)
            elif parameter == "plot_FF_mag":
                valid = True ; plot_FF_mag = str2bool(value)
            elif parameter == "plot_FF_z":
                valid = True ; plot_FF_z = str2bool(value)
            elif parameter == "plot_FF_x":
                valid = True ; plot_FF_x = str2bool(value)
            elif parameter == "plot_FF_y":
                valid = True ; plot_FF_y = str2bool(value)

#-------------------------------------------------------------------------------
#           Output velocity vector and plot scale specifications
#-------------------------------------------------------------------------------

            elif parameter == "vector_overlay":
                valid = True ; vector_overlay = str2bool(value)
            elif parameter == "vector_color":
                valid = True ; vector_color = str(value)
            elif parameter == "vector_length":
                valid = True ; vector_length = float(value)
            elif parameter == "vector_width":
                valid = True ; vector_width = float(value)
            elif parameter == "vector_interval":
                valid = True ; vector_interval = int(value)
            elif parameter == "annulus_interval":
                valid = True ; annulus_interval = int(value)
            elif parameter == "clip_rho":
                valid = True ; clip_rho = float(value)
            elif parameter == "gpow_rho":
                valid = True ; gpow_rho = float(value)
            elif parameter == "clip_T":
                valid = True ; clip_T = float(value)
            elif parameter == "gpow_T":
                valid = True ; gpow_T = float(value)
            elif parameter == "fix_scale_T":
                valid = True ; fix_scale_T = str2bool(value)
            elif parameter == "clip_u":
                valid = True ; clip_u = float(value)
            elif parameter == "gpow_u":
                valid = True ; gpow_u = float(value)
            elif parameter == "clip_nu_total":
                valid = True ; clip_nu_total = float(value)
            elif parameter == "gpow_nu_total":
                valid = True ; gpow_nu_total = float(value)
            elif parameter == "clip_Tau":
                valid = True ; clip_Tau = float(value)
            elif parameter == "gpow_Tau":
                valid = True ; gpow_Tau = float(value)

#-------------------------------------------------------------------------------
#           Performance specifications
#-------------------------------------------------------------------------------

            if parameter == "code":
                valid = True ; code = str(value)
            elif parameter == "overlap":
                valid = True ; overlap = str2bool(value)
            elif parameter == "nxi_min":
                valid = True ; nxi_min = int(value)
            elif parameter == "nxi_max":
                valid = True ; nxi_max = int(value)
            elif parameter == "nyi_min":
                valid = True
                if value != "unspecified": nyi_min = int(value)
                else:                      nyi_min = str(value)
            elif parameter == "nyi_max":
                valid = True
                if value != "unspecified": nyi_max = int(value)
                else:                      nyi_max = str(value)
            elif parameter == "nzi_min":
                valid = True
                if value != "unspecified": nzi_min = int(value)
                else:                      nzi_min = str(value)
            elif parameter == "nzi_max":
                valid = True
                if value != "unspecified": nzi_max = int(value)
                else:                      nzi_max = str(value)
            elif parameter == "domain_decomposition":
                valid = True ; domain_decomposition = str(value)

#-------------------------------------------------------------------------------
#           Print a warning any lines which do not set valid paramaters
#-------------------------------------------------------------------------------

            if not valid:
                nwarning += 1
                if rank == 0:
                    for iout in range(1):
                        if   iout == 0: output_print = sys.stderr
                        elif iout == 1: output_print = sys.stdout
                        print('WARNING: line number = %d, line = "%s"' %(iline,line.rstrip()),file=output_print,flush=True)
                        print("         %s is an invalid parameter name, this line is being ignored" %(parameter),file=output_print,flush=True)

    if nerror > 0:
        if rank == 0:
            for iout in range(2):
                if   iout == 0: output_print = sys.stderr
                elif iout == 1: output_print = sys.stdout
                print(file=output_print)
                print("Errors were detected in the input file which need to be fixed to do a run",file=output_print,flush=True)
                print(file=output_print,flush=True)
                print("       Valid input lines must be of form: parameter_name = parameter_value",file=output_print,flush=True)
                print("       where starting from #, input is considered a comment and is ignored",file=output_print,flush=True)
                print("",file=output_print)
                print("Exited run!",file=output_print)
                print("==================================================================",file=output_print,flush=True)
        exit()

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

