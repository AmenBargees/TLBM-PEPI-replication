
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

def str2bool(value,fouts):

    global parameter, iline, line, n_warnings

    if   value == "True"  or value == "true":
        result = True
    elif value == "False" or value == "false":
        result = False
    else:
        result = False
        n_warnings += 1
        if rank == 0:
            for i in range(len(fouts)):
                print('WARNING: line number %d, line: "%s"' %(iline,line.rstrip())                ,file=fouts[i],flush=True)
                print('         Must specify True or False for boolean variable "%s"' %(parameter),file=fouts[i],flush=True)
                print('         "%s" is being set to False' %(parameter)                          ,file=fouts[i],flush=True)
    return result

def get_parameters(filename,fouts=[],rank=0):

    global parameter, iline, line, n_warnings

#   Test option

    global test_option

#   Run title and description for plots

    global run_title, description, description_color, description_box, description_location, description_alignment, description_text_x

#   Model size, time steps and geometry

    global L_p, nx, ny, nz, nt, geometry, theta1, theta2

#   Physical specifications

    global state, pturb, Ra, Pr, R_p, L_p, rho_p, eta_p, ka_p, g_p, beta_p, T0_p, dT_p
    global rho_z_file, P_z_file, T_z_file
    global T_p_random, dT_p_random, ncells_x, ncells_y, cell_amplitude_body, cell_amplitude_boundaries, cell_power, amp_sinusoid, x_periods, y_periods, z_periods
    global Tdep_viscosity, gamma_T, E_offset, rheology, Tau_y0_p, Tau_yi_p, N, TTG_partial_melting

#   Boundary condition specifications

    global BC_sides, BC_upper, BC_lower

#   Physical time step parameters
#       If g  is specified, nu will be calculated from g
#       If nu is specified, g  will be calculated from nu
#       If both nu and g are specified, nu will take precedence

    global dt_p, nu_l, optimize_nu, g_l, Ma_sim, u_p_geo

#   Mach control and accuracy/stability parameters

    global max_mach_limit, min_mach_limit, Inc_vary_dt, start_vary_dt
    global factor_increase_dt,factor_decrease_dt,factor_max_vary_dt,mach_limit_to_exit
    global nu_l_min

#   Accuracy parameters

    global float_bits, edge_bits, data_bits, shifted_distribution

#   Output specifiations

    global output_plots, output_data, output_bmk, output_units

#   Output increment parameters

    global Inc_mach, Inc_plots, Inc_data, Inc_bmk, Inc_performance

#   Output plot parameters

    global plane, coordinates
    global plot_size, plot_size_factor, plot_colorbar, colormap, latex
    global plot_rho, plot_T, plot_u, plot_vectors, plot_nu_total, plot_Tau
    global plot_rho_eff, plot_rho_s, plot_M, plot_ce_pd, plot_T_sol, plot_T_liq, plot_P

#   Output velocity vector and plot scale parameters

    global vector_overlay, vector_color, vector_length, vector_width, vector_interval, annulus_interval
    global clip_rho, gpow_rho, clip_T, gpow_T, fix_scale_T, clip_u, gpow_u
    global clip_nu_total, gpow_nu_total, clip_Tau, gpow_Tau

#   Performance parameters

    global code, overlap, domain_decomposition
    global nxi_min, nxi_max, nyi_min, nyi_max, nzi_min, nzi_max
    global run_type

#-------------------------------------------------------------------------------
#   Get read and set parameters
#-------------------------------------------------------------------------------

    if rank == 0:
        for i in range(len(fouts)):
            print("Reading input file ...",end="",file=fouts[i],flush=True)

    file_input  = open(filename,"r")
    
    iline      = 0       # Line number being read
    nerror     = 0       # Number of input errors
    n_warnings = 0       # Number of input warnings
    verbose    = False   # Verbose mode prints out each parameter as it is set

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
                    for i in range(len(fouts)):
                        if nerror == 1: print()
                        print('ERROR: line number %d, line: "%s"' %(iline,line.rstrip()),file=fouts[i],flush=True)
        elif len(parameters) == 3:
            parameter = parameters[0]
            value     = parameters[2]
            valid     = False
            if rank == 0:
                for i in range(len(fouts)):
                    if verbose:
                        if nerror == 1: print()
                        print("Setting parameter ",parameter," to ",value               ,file=fouts[i],flush=True)

#-------------------------------------------------------------------------------
#           Test option
#-------------------------------------------------------------------------------

            if parameter == "test_option":
                valid = True ; test_option = value

#-------------------------------------------------------------------------------
#           Run title and description for plots
#-------------------------------------------------------------------------------

            elif parameter == "run_title":
                valid = True ; run_title = value
            elif parameter == "description":
                valid = True ; description = value
            elif parameter == "description_color":
                valid = True ; description_color = str(value)
            elif parameter == "description_box":
                valid = True ; description_box = str2bool(value,fouts)
            elif parameter == "description_location":
                valid = True ; description_location = str(value)
            elif parameter == "description_alignment":
                valid = True ;
                if   value == "left":   description_alignment = "left"
                elif value == "center": description_alignment = "center"
                elif value == "right":  description_alignment = "right"
                else:                   description_alignment = "center"

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
                valid = True
                if value!= "unspecified": nz = int(value)
                else:                     nz = str(value)
            elif parameter == "nt":
                valid = True ; nt = str(value)
            elif parameter == "geometry":
                valid = True ; geometry = str(value)
            elif parameter == "theta1":
                valid = True ; theta1 = float(value)
            elif parameter == "theta2":
                valid = True ; theta2 = float(value)

#-------------------------------------------------------------------------------
#           Physical specifications
#-------------------------------------------------------------------------------

            elif parameter == "state":
                valid = True
                state = str(value)
            elif parameter == "pturb":
                valid = True
                if type(value) == str: pturb = str(value)
                else:                  pturb = float(value)
            elif parameter == "Ra":
                valid = True
                if value != "unspecified": Ra = float(value)
                else:                      Ra = str(value)
            elif parameter == "Pr":
                valid = True
                if value != "unspecified": Pr = float(value)
                else:                      Pr = str(value)
            elif parameter == "rho_p":
                valid = True ; rho_p = float(value)
            elif parameter == "rho_z_file":
                valid = True ; rho_z_file = str(value)
            elif parameter == "P_z_file":
                valid = True ; P_z_file = str(value)
            elif parameter == "T_z_file":
                valid = True ; T_z_file = str(value)
            elif parameter == "eta_p":
                valid = True
                if value != "unspecified": eta_p = float(value)
                else                     : eta_p = value
            elif parameter == "ka_p":
                valid = True
                if value != "unspecified": ka_p = float(value)
                else                     : ka_p = str(value)
            elif parameter == "g_p":
                valid = True ; g_p = float(value)
            elif parameter == "beta_p":
                valid = True ; beta_p = float(value)
            elif parameter == "T0_p":
                valid = True ; T0_p = float(value)
            elif parameter == "dT_p":
                valid = True ; dT_p = float(value)
            elif parameter == "T_p_random":
                valid = True ; T_p_random = float(value)
            elif parameter == "dT_p_random":
                valid = True ; dT_p_random = float(value)
            elif parameter == "ncells_x":
                valid = True ; ncells_x = int(value)
            elif parameter == "ncells_y":
                valid = True ; ncells_y = int(value)
            elif parameter == "cell_amplitude_boundaries":
                valid = True ; cell_amplitude_boundaries = float(value)
            elif parameter == "cell_amplitude_body":
                valid = True ; cell_amplitude_body = float(value)
            elif parameter == "cell_power":
                valid = True ; cell_power = float(value)
            elif parameter == "amp_sinusoid":
                valid = True ; amp_sinusoid = float(value)
            elif parameter == "x_periods":
                valid = True ; x_periods = int(value)
            elif parameter == "y_periods":
                valid = True ; y_periods = int(value)
            elif parameter == "z_periods":
                valid = True ; z_periods = int(value)
            elif parameter == "Tdep_viscosity":
                valid = True ; Tdep_viscosity = str(value)
            elif parameter == "gamma_T":
                valid = True ; gamma_T = float(value)
            elif parameter == "E_offset":
                valid = True ; E_offset = float(value)
            elif parameter == "rheology":
                valid = True ; rheology = str(value)
            elif parameter == "Tau_y0_p":
                valid = True ; Tau_y0_p = float(value)
            elif parameter == "Tau_yi_p":
                valid = True ; Tau_yi_p = float(value)
            elif parameter == "N":
                valid = True ; N = float(value)
            elif parameter == "TTG_partial_melting":
                valid = True ; TTG_partial_melting = str(value)

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
#               if dt_p is specified, it takes overall precedence
#
#-------------------------------------------------------------------------------

            elif parameter == "dt_p" or parameter == "time_per_timestep": # Retain back compatibility of time_per_timestep
                valid = True
                dt_p = str(value)
            elif parameter == "Ma_sim":
                valid = True
                if value!= "default" and value != "unspecified": Ma_sim = float(value)
                else:                                            Ma_sim = str(value)
            elif parameter == "u_p_geo":
                valid = True
                if value!= "default" and value != "unspecified": u_p_geo = float(value)
                else:                                            u_p_geo = str(value)
            elif parameter == "nu_l":
                valid = True
                if value != "unspecified": nu_l = float(value)
                else:                      nu_l = str(value)
            elif parameter == "optimize_nu":
                if value == "True" or value == "False":
                    valid = True ; optimize_nu = str2bool(value,fouts)
                else:
                    valid = True ; optimize_nu = float(value)
            elif parameter == "g_l":
                valid = True
                if value != "unspecified": g_l = float(value)
                else:                      g_l = str(value)

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
                valid = True ; shifted_distribution = str2bool(value,fouts) ;

#-------------------------------------------------------------------------------
#           Output specifiations
#-------------------------------------------------------------------------------

            elif parameter == "output_plots":
                valid = True ; output_plots  = str2bool(value,fouts)
            elif parameter == "output_data":
                valid = True ; output_data   = str2bool(value,fouts)
            elif parameter == "output_bmk":
                valid = True ; output_bmk    = str2bool(value,fouts)
            elif parameter == "output_units":
                valid = True ; output_units  = str(value)

#-------------------------------------------------------------------------------
#           Output increment specifications
#-------------------------------------------------------------------------------

            elif parameter == "Inc_mach":
                valid = True ; Inc_mach        = str(value)
            elif parameter == "Inc_plots":
                valid = True ; Inc_plots       = str(value)
            elif parameter == "Inc_data":
                valid = True ; Inc_data        = str(value)
            elif parameter == "Inc_bmk":
                valid = True ; Inc_bmk         = str(value)
            elif parameter == "Inc_performance":
                valid = True ; Inc_performance = str(value)

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
                valid = True ; plot_colorbar = str2bool(value,fouts)
            elif parameter == "colormap":
                valid = True ; colormap = str(value)
            elif parameter == "latex":
                valid = True ; latex = str2bool(value,fouts)
            elif parameter == "plot_rho":
                valid = True ; plot_rho = str2bool(value,fouts)
            elif parameter == "plot_rho_eff":
                valid = True ; plot_rho_eff = str2bool(value,fouts)
            elif parameter == "plot_rho_s":
                valid = True ; plot_rho_s = str2bool(value,fouts)
            elif parameter == "plot_M": 
                valid = True ; plot_M = str2bool(value,fouts)
            elif parameter == "plot_ce_pd":
                valid = True ; plot_ce_pd = str2bool(value,fouts)
            elif parameter == "plot_T_sol":
                valid = True ; plot_T_sol = str2bool(value,fouts)
            elif parameter == "plot_T_liq":
                valid = True ; plot_T_liq = str2bool(value,fouts)
            elif parameter == "plot_P":
                valid = True ; plot_P = str2bool(value,fouts)
            elif parameter == "plot_T":
                valid = True ; plot_T = str2bool(value,fouts)
            elif parameter == "plot_u":
                valid = True ; plot_u = str2bool(value,fouts)
            elif parameter == "plot_vectors":
                valid = True ; plot_vectors = str2bool(value,fouts)
            elif parameter == "plot_nu_total":
                valid = True ; plot_nu_total = str2bool(value,fouts)
            elif parameter == "plot_Tau":
                valid = True ; plot_Tau = str2bool(value,fouts)

#-------------------------------------------------------------------------------
#           Output velocity vector and plot scale specifications
#-------------------------------------------------------------------------------

            elif parameter == "vector_overlay":
                valid = True ; vector_overlay = str2bool(value,fouts)
            elif parameter == "vector_color":
                valid = True ; vector_color = str(value)
            elif parameter == "vector_length":
                valid = True ; vector_length = float(value)
            elif parameter == "vector_width":
                valid = True ; vector_width = float(value)
            elif parameter == "vector_interval":
                valid = True
                if value == "default": vector_interval  = str(value)
                else                 : vector_interval  = int(value)
            elif parameter == "annulus_interval":
                valid = True
                if value == "default": annulus_interval = str(value)
                else                 : annulus_interval = int(value)
            elif parameter == "clip_rho":
                valid = True ; clip_rho = float(value)
            elif parameter == "gpow_rho":
                valid = True ; gpow_rho = float(value)
            elif parameter == "clip_T":
                valid = True ; clip_T = float(value)
            elif parameter == "gpow_T":
                valid = True ; gpow_T = float(value)
            elif parameter == "fix_scale_T":
                valid = True ; fix_scale_T = str2bool(value,fouts)
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
                valid = True ; overlap = str2bool(value,fouts)
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
            elif parameter == "run_type":
                valid = True ; run_type             = str(value)

#-------------------------------------------------------------------------------
#           Print a warning any lines which do not set valid paramaters
#-------------------------------------------------------------------------------

            if not valid:
                n_warnings += 1
                if rank == 0:
                    for i in range(len(fouts)):
                        print('WARNING: line number = %d, line = "%s"' %(iline,line.rstrip())                    ,file=fouts[i],flush=True)
                        print("         %s is an invalid parameter name, this line is being ignored" %(parameter),file=fouts[i],flush=True)

    if nerror > 0:
        if rank == 0:
            for i in range(len(fouts)):
                print(file=fouts[i])
                print("Errors were detected in the input file which need to be fixed to do a run"                 ,file=fouts[i],flush=True)
                print(                                                                                             file=fouts[i],flush=True)
                print("       Valid input lines must be of form: parameter_name = parameter_value"                ,file=fouts[i],flush=True)
                print("       where starting from #, input is considered a comment and is ignored"                ,file=fouts[i],flush=True)
                print(                                                                                             file=fouts[i],flush=True)
                print("Exited run!"                                                                               ,file=fouts[i],flush=True)
                print("=================================================================="                        ,file=fouts[i],flush=True)
        exit()

#-------------------------------------------------------------------------------
#   Print the status if the file was successfully read
#-------------------------------------------------------------------------------

    if rank == 0:
        if n_warnings > 0:
            for i in range(len(fouts)):
                print(" file read with %d warnings." %(n_warnings)                                                ,file=fouts[i],flush=True)
        else:
            for i in range(len(fouts)):
                print(" read successfully!"                                                                       ,file=fouts[i],flush=True)

