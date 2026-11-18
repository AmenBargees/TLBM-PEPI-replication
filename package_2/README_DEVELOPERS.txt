============================================================================================
DEVELOPERS
============================================================================================
This document is a brief overview for developers describing how to add new input parameters,
new physics, and new plot & data output options to GTXLBM.

IMPORTANT: Changes should be made to GTXLBM.py, not GTXLBM_run.py which is a copy of GTXLBM.py
used for runs. When you are ready to run and test changes, GTXLBM_run.py can be updated via
the update script which overwrites the current version of GTXLBM_run.py and inserts information
including the GTXLBM version and its release date, your user name, and the date of this update.
--------------------------------------------------------------------------------------------
(1) New input parameters
--------------------------------------------------------------------------------------------

    (a) New input parameters require making additions of the new parameters to the get_parameters() function
        of the get_parameters.py file. There are two changes required: (i) Add the new parameter
        to the global variables of get_parameters(), (ii) Add lines to read the new parameter.
        With these changes, both GTXLBM.py and plot_npy.py can read the new parameter from
        the input file such as the line below to set the value of a new parameter called "new_var".

        new_var  = 2.34e5  # Be sure to add a description of new_var in your input file & its units in MKS

       For example, to add a new float parameter called "new_var", the following changes are required.

       (i)  Add "new_var" to the global parameters list in a logical location. For example, if "new_var"
            is a physical specification, add new_var to the globals as below:

            #   Physical specifications

                global Ra_p, ....
                ...
                global new_var

       (ii) Add lines in the logical location to read the new input parameter. For example, if "new_var"
            is a physical specification and a float, add "new_var" as below:

            #-------------------------------------------------------------------------------
            #  Physical specifications
            #-------------------------------------------------------------------------------

            elif parameter == "Ra_p":
                valid = True
                if value != "unspecified": Ra_p = float(value)
                else:                      Ra_p = str(value)
            ...
            elif parameter == "new_var":
                valid = True ; new_var = float(value)

            NOTE: If "new_var" is a bool, string or int, follow the same protocol as above, but use
            new_var = str2bool(value), new_var = str(value), or new_var = int(value) as appropriate
            instead of new_var = float(value)

    (b) New input parameters also require changes to the master file TLBM_master.py, followed by running
        the update.s script to transcribe these changes to the run code, GTXLBM.py. There are a number
        of changes required which are: (i) Initialize the new parameter to some reasonable default, (ii)
        Print the value of the new parameter that was read from the input file in TLBM_master.py. For
        example, consider the example of the new parameter "new_var" which is a float and a physical
        specification. Changes are required as shown below:

        (i)  To add new parameter "new_var", a physical specification of type float, add new_var to the
             lines as below.

             #-------------------------------------------------------------------------------
             # Physical specifications
             #-------------------------------------------------------------------------------

             Ra_p                 = 1.e5          # Physical Rayleigh number
             ...
             eta_p                = 5e21          # Dynamic viscosity               (Pa.s = kg/(m.s))
             ...
             N                    = 4.0           # Exponent for plastic rheology equation: Tau = (Tau_y/Tau)^N
             new_var              = 1.23e4        # Put a description of the new parameter here & specify its MKS units

        (ii) To print the new parameter "new_var", a physical specification of type float, in the output file,
             add new_var as below:

             print("------------------------------------------------------------------",flush=True)
             print("Physical specifications",flush=True)
             print("------------------------------------------------------------------",flush=True)
             print("Ra_p                         = %8.3g" %(Ra_p),flush=True)
             ...
             print("eta_p       (kg/(m.s)        = %8.3g" %(eta_p),flush=True)
             print("N           (exponent)       = ",N,end="",flush=True)
             print("new_var     (MKS units here) = %8.3g" %(new_var),flush=True)

--------------------------------------------------------------------------------------------
(3a) New physics, step 1: Add new physics to Python code and verify it!
--------------------------------------------------------------------------------------------

If new physics is required, the best way to do this is to first follow steps (1) and (2) above to
get any new parameters that are required and ensure these are printed to the output files to record
and verify what a give run is doing, and then to make changes to the Python code first in the
appropriate place in the code. Once the Python code has been verified to do what is desired such
as by analyzing output data and plots, the new physics can be added to the C code to obtain higher
performance, typically at least 3x and up to around 10x. Below, is an artificial example where
a new temperature dependent viscosity is added into the Python code to illustrate. Suppose that instead
of Reynolds or Arrhenius temperature dependence which are already implemented (ie. temperature = reynolds,
gamma_T = value_for_power or rheology = archean2, temperature = arrhenius), that we want a power law
temperature dependence of form: nu = nu_0 (T/T_ref)^gamma2_T. In this case, add gamma2_T and T_ref to
the input parameters as shown above (see example of new_var), and make additions as below to
the Python code. Note: "temperature = power_law" and "gamma2_T = -2.34" must be set in your input
file assuming you want to invoke the new power law with gamma2_T = -2.34 per the code below.

#           Calculate the linear viscosity term nu_linear

            if temperature == "arrhenius":                  # Arrhenius temperature dependence
                
                act_E = 270e3   # Activation energy         units: J/mol
                act_V = 5e-6    # Activation volume         units: m3/mol
                gas_C = 8.314   # Gas constant              units: J/(mol*K)
                ref_K = 1600    # Reference temperature     units: K

                arg = (act_E + P*1e9 * act_V)/(gas_C * K) - (act_E + P*1e9 * act_V)/(gas_C * ref_K)
                arg = np.where(arg<87, arg, 87)
                nu_linear = nu_f * np.exp(arg)

            elif gamma_T > 0 and temperature == "reynolds": # Reynolds temperature dependence

                nu_linear = nu_f * np.exp(-gamma_T*eta)

            elif gamma2_T != 0 and temperature == "power_law": # Power law temperature dependence

                T         = eta/DR_2 + T0_p + dT_p/2           # Calculate temperature from energy density
                nu_linear = nu_f * np.power(T/T_ref,gamma2_T)  # New power law temperature dependence

            else:                                              # No temperature dependence

                nu_linear = nu_f

--------------------------------------------------------------------------------------------
(3a) Adding new physics, step 2: Add new physics to C code and verify it gives the same
                                 result as the Python code.
--------------------------------------------------------------------------------------------
Adding the calculations for new physics requires the equivalent changes to the C code as
were made to the Python code if you want to get the large performance gains of using
optimized C code with cache blocking, typically 3x and up to around 10x performance gains.
This involves: (a) specify any new input arguments to the C code in the Python argtypes
statement in TLBM_master.py as appropriate, (b) specify any new input arguments to the
C function in GTXLBM.c, (c) ensure any passed variables or arrays to the C function are
of type np_float for floats and np_int for ints, (d) add the new physics to the C kernel
routine GTXLBM_kernel.c. As an example, the power law temperature dependence as above
can be added by the following three steps:

(a) Change the Python code as below to specify the five new required float arguments to the GTXLBM C function.
    It is important to make these changes as below. Namely, in two places which are: (1) to the commented section
    before the argtypes specification of GTXLBM C function so you can see the names of the arguments to
    the GTXLBM() C function and exactly where they are in the argument list, and (2) the argtypes
    specification line GTXLBM.py where the new argument types must be specified. Namely, five new
    arguments - "gamma2_T, T_ref, DR_2, T0_p & dT_p" - were added before the last argument named "rank",
    as these will be added to the arguments of the C code in the same location, ie. before "rank".
    See "C_FLOAT, C_FLOAT, C_FLOAT, C_FLOAT, C_FLOAT" before the last argument of type "c_int" which
    corresponds to the last argument named "rank".

    # The arguments below show the meaning of the arguments in the argtypes specifications
    #
    #  GTXLBM(
    #  na,D,nx_i,ny_i,nz_i,mx_i,my_i,mz_i,xedge,yedge,zedge, \
    #  ...
    #  gamma2_T,T_ref,DR_2,T0_p,dT_p, \
    #  rank
    #  )

    GTXLBM_argtypes = [ \
        c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int, c_int, \
        ...
        C_FLOAT, C_FLOAT, C_FLOAT, C_FLOAT, C_FLOAT, \
        c_int \
    ]

(b) Specify the five new input arguments to the C function in GTXLBM.c in the same location
    as in the argtypes statement in section (1) as illustrated below. This needs to be done in two locations,
    once for the 32 bit function and the same change for the 64 bit function.

    ...
    void GTXLBM_32(
        int na, int D, int nx, int ny, int nz, int mx, int my, int mz, int xedge, int yedge, int zedge,
        ...
        FLOAT gamma2_T, FLOAT T_ref, FLOAT DR_2, FLOAT T0_p, FLOAT dT_p,
        int rank
    )
    ...
    void GTXLBM_64(
        int na, int D, int nx, int ny, int nz, int mx, int my, int mz, int xedge, int yedge, int zedge,
        ...
        FLOAT gamma2_T, FLOAT T_ref, FLOAT DR_2, FLOAT T0_p, FLOAT dT_p,
        int rank
    )

    Note: FLOAT is defined as "float" for the 32 bit C function gTLMBX(), and as "double"
    for the 64 bit C function gTLMBX_64().

(c) Ensure passed variables or arrays to the C function are of type np_float for floats or np_int for ints
    where np_float = np.float32 or np.float64 depending on input parameter "float_bits", and np_int = np.int32.
    If you need to pass a new string to the C function, it must be converted to the C "char" type as shown below
    the Python string named "temperature" is converted to the Python name "Temperature" (bytes) which is of the
    same type as the C "char" type. Hence, "Temperature" can be passed to C, but "temperature" cannot.

    ...

    # Initialize float and int constants to be of the correct type

    beta_g      = np_float(-3 * beta_p * g_l)        # Boussinesq factor
    ...
    gamma2_T = np_float(gamma2_T)                    # Power for power law temperature dependence
    T_ref    = np_float(T_ref)                       # Reference temperature, power law temperature dependence
    DR_2     = np_float(DR_2)                        # Conversion factor energy density to temperature
    T0_p     = np_float(T0_p)                        # Surface temperature
    dT_p     = np_float(dT_p)                        # Temperature difference between upper & lower boundaries
    ...

    # Convert Python strings to bytes which is equivalent to the C "char" type which
    # can be passed to the C GTXLBM() function

    if   temperature == "reynolds":  Temperature = b"reynolds"

    ...

(d) Add the new physics to the C kernel routine GTXLBM_kernel.c as illustrated below for the new
    power law temperature dependence example above.

    ...
    FLOAT nu_linear[nx], nu_plastic[nx], nu_total[nx], tau_f[nx], gi[nx];
    FLOAT T[nx]; // Declare a new 1D temperature array over the cache blocking dimension = last dimension = x-axis
    ...
    //---------------------------------------------------------------
    //              Macroscopic properties ...
    //---------------------------------------------------------------
    ...
    for(x=startx;x<endx;x+=incx) {
        eta[y][z][x] /= rho [y][z][x];
        T        [x]  = (eta[y][z][x]/DR_2+T0_p+dT_p/2); // Calculate the temperature from energy density
    }
    ...
    //---------------------------------------------------------------
    //  Calculate nu_linear, so any temperature dependence
    //---------------------------------------------------------------

        if(strcmp(temperature,"reynolds") == 0 && gamma_T > 1.e-10)  {      # Reynolds      temperature dependence
            ....
        }
        ....
        else if(strcmp(temperature,"power_law") == 0 && gamma2_T != 0.0)  { # New power law temperature dependence
            for(x=startx;x<endx;x+=incx) {
                nu_linear[x] = nu_f[z] * POW(-gamma2_T*T[x]/T_ref);
            }
        }
        else {                                                              # No            temperature dependence
            ....
        }

--------------------------------------------------------------------------------------------
(4) New output plot and data file output options
--------------------------------------------------------------------------------------------

Any new variables you wish to save data files and/or plot require TLBM_master.py and plot_npy.py
to be updated to write/read data files of the new variables, and to plot the new variables. There
are three additions to make: (a) Add new bool variable names to allow the new variable to be specified
for plotting, (b) Add lines to TLBM_master.py to save data files and plot new variables, and
(c) Add lines to plot_npy.py to allow it to read and plot data files for new variables.

(a) Follow the same approach as above explaining how to add new parameters, but this time, you
    will need to add a bool variable such as "plot_new_var", and if you want it to have its own
    plot scale parameters (clip and gpow), then you will also need to add "clip_new_var" and
    gpow_new_var.

(b) Add lines of code into TLBM_master.py at the appropriate location as below to the save data files
    new_var_*.npy and plot the new variable to new_var_*.png in the data and plot directories
    for the project and run.

    # Write and plot temperature

      write_and_plot_field(eta/DR_2,"T",r'$\displaystyle T$',"Degrees (C)",r'$^\circ C$', \
          save=True,plot=plot_T,clip=clip_T,gpow=gpow_T,vmin=T0_p,vmax=T0_p+dT_p,c0=T0_p+dT_p/2,c1=1)
    ...

    # Write and plot variable new_var

    if plot_new_var:
        write_and_plot_field(new_var,"new_var",r'$\displaystyle new_{var}$',"Units_new_var",r'Units_new_var',save=True,plot=plot_new_var)
    ...

    NOTE: Required arguments of write_and_plot_field() are:
            1st argument = variable name on cores (must be a of type np_float)
            2nd argument = string to identify new_var files and the plain text label for new_var on plots,
            3rd argument = latex quality label on plots for new_var (only used for if "latex = True" is specified),
            4th argument = plain text string for labels on the colorbar to specify the units of new_var
            5th argument = latex quality label on the colorbar to specify the units of new_var (only for "latex = True)
          Optional arguments are:
            save           True or False to specify whether to save an npy data file, plot = True or False to specify
                           whether to save a .png plot file
            vmin & vmax    The min and max values to plot if specified. If vmin and vma are not specified,
                           the plot will show the full range from the minimum to the maximum of the variable
          c0 and c1        Constants to use to rescale the values of new_var according to the formula:
                           new_var_for_plot = c0 + c1 * new_var. This allows certain lattice values to be
                           rescaled to physical values such as converting energy density eta to temperature.
                           If unspecified, c0 defaults to 0 and c1 to 1 which means no rescaling.

(c) Add lines to plot_npy.py to allow it to read and plot data files for new variables as below.

    #-------------------------------------------------------------------------------
    # Plot temperature
    #-------------------------------------------------------------------------------

    if plot_T:
        read_and_plot_field("T",r'$\displaystyle T$',"Degrees (C)",r'$^\circ C$', \
                clip=clip_T,gpow=gpow_T,vmin=T0_p,vmax=T0_p+dT_p,c0=T0_p+dT_p/2,c1=1)
    ...

    #-------------------------------------------------------------------------------
    # Plot new_var
    #-------------------------------------------------------------------------------

    if plot_new_var:
        read_and_plot_field("new_var",r'$\displaystyle new_{var}$',"Units of new_var",r'Units of new_var')

    NOTE: Arguments of read_and_plot_field() are the same as write_and_plot_field() but without the 1st argument
    (ie. without a variable name) as this function reads data files rather than getting the field from the domains
    of the variable on each core.

--------------------------------------------------------------------------------------------
(5) Initializing temperature or pressure profiles, or other variables
--------------------------------------------------------------------------------------------

In some cases, it may be necessary to initialize a space or depth dependent variable. This
can be done by adding an input parameter with a filename, and then reading the variable from
the file, and then putting this space dependent variable onto the cores. For example,
to initialize a field P0 on the cores (ie. a field or array on the domains) corresponding to
a depth dependent pressure, add a new input parameter "P_z_file" which specifies the filename
of numpy (ie. .npy) file that contains the pressure profile. The file can be created
using numpy.save() of a 1D array with the same z-dimension as specified in the input file.
The file must be created and put into ./projects/project/run/, and it can be read within
GTXLBM.py and put onto the cores as below. This is how rheology == "archean2" gets
the temperature, density and pressure profiles.

     P_z_file   = run_dir + P_z_file   + ".npy"
     P_z        = np_float(np.load(P_z_file))
     P0         = put_field_z_2_xyz(P_z,default=0)
     del P_z

If you wish to initialize a field that is a function of x, y and z, the procedure is the same
as above but instead of using put_field_z_2_xyz() - which puts a 1D function of z onto a 3D
field on the cores - use use put_field() which puts a 3D field onto a 3D field the cores.
And of course, the numpy file being read must have the same dimensions as the 3D field,
so of dimensions such as via being initialized by np.zeros((ny,nz,nx),dtype=np.float32).

Note: The value of the default in the put_field() and put_field_z_2_xyz() functions is put
onto the edges of domains and has no effect on the TLBM calculations themselves. This is
because edges are extra space in the arrays on the cores (domains) to allow enough space
for the getting of edges from adjacent domains which is required for the streaming step
of the LBM.  However, in Python, as calculations are done on the entire field including
edges, it is useful to initialize the edges to a suitable value as this can avoid overflows
in Python for the calculations on the edges, even if such overflows can be ignored as
they are irrelevent to the TLBM calculations.

--------------------------------------------------------------------------------------------
Calculations in Python or C
--------------------------------------------------------------------------------------------

GTXLBM.py has an option to do the calculations in Python or C. As such, new physics
should initially be implemented and tested in Python which easy, and when the Python
implementation is verified to be correct, it is relatively easy to make the same changes
to the C code kernel code GTXLBM_kernel.c which uses cache blocking and is highly
optimized. This can lead to a three to ten fold increase in the performance, so it is
well worthwhile, especially when doing larger high resolution 2D runs or 3D runs.

--------------------------------------------------------------------------------------------
Accuracy and stability
--------------------------------------------------------------------------------------------

The TLBM involves relaxing the number densities carrying mass to the Boltzmann equilibrium
which is a Taylors expansion to second order in fluid velocity. The terms in the expansion
have (u/c_s)^n where n is 1 and 2, and c_s is the sound speed in the lattice, so powers of
the mach number. Hence, the TLBM requires the mach number to be small enough relative to
achieve a high accuracy. Furthermore, if the mach number becomes too high, the TLBM can
blow up (become unstable) leading to to infinities and nans. Based on modified Blankenbach
benchmarks for finite & infinite Pr with PDE solvers [4], the TLBM solution is highly
accurate for mach number <~ 0.01, and TLBM solutions with Pr=1000 solutions precisely match
with Pr=infinity PDE solutions at such low mach numbers (ie. the present day earth with
Pr >~ 10^30 ~ infinity can be modeled with with the finite Pr number TLBM using Pr ~ 1000).
In some cases, 64 bit precision may be required such as for high Pr runs which may
lead to very low lattice thermal diffusivity, and hence, relaxation times for distribution
g carrying energy density very close to 0.5. ie. tau_g = 0.5 + 3 nu_g = 0.5 + 3 kappa_l
where kappa_l is the lattice thermal diffusivity.

Note that one can consider that the mach number plays a similar role as the Courant condition
in numerical solutions to PDEs. For example, the Courant condition may be that the timestep
dt <= s.dx/v for stability where v is the maximum velocity (eg. P wavespeed), and s is a
number less than unity such as s = 0.5, but s << 1 is generally required for accuracy (eg.
s = 0.1, so dt ~ 0.1 dx/v).

Similarly, u = flow velocity in the TLBM must be significantly less than the mach number, so
u/c_s must be small for stability, and smaller still for accuracy such as mach number u/c_s
<~ 0.01 or <~ 0.1 depending on the desired accuracy. As the speed of sound in the LBM is c_s =
1/sqrt(3), and dx = dt = 1 (in the standard LBM), this means sqrt(3).u << dx/dt or dt << dx/(sqrt(3).u)
or dt < s_LBM dx/u where s_LBM = 0.01 (highly accurate), s_LBM = 0.04 (less accurate but
probably Ok for various problems), and s_LBM = 0.2 (not so accurate and the LBM will model
partially compressible flow, but may be Ok for some problems). Hence, the mach condition for TLBM
runs has the same form as the Courant condition for PDE solvers.

In geodynamics problems or CFD in general, the fluid flow velocity may vary enormously for a
given set of simulation parameters, and may change during a simulation.

As such, GTXLBM.py has variable timestepping to adapt to the flow velocity to ensure
a desired accuracy is maintained, while minimizing the runtime, see below! This technique
increases the timestep when the flow is slow, and decreases the timestep when the flow
is fast. This maintains precision and minimizes computational time, no matter the variations
in the flow speed throughout the simulation, or what initial parameter choices were made
which control the initial physical time step (see reference [3] for the time scaling
formula which relates LBM parameters to the physical timestep).

--------------------------------------------------------------------------------------------
Minimizing run time
--------------------------------------------------------------------------------------------

GTXLBM.py has variable time stepping which adjusts the LBM time step to keep the mach
number within a specified range. While benchmarks [4] indicate that very accurate solutions
can be achieved for mach number ~ 0.01, subsequent tests reveal that runs at at much higher
mach numbers may be sufficiently accurate (ie. the qualitative behavior and phenomena
one observes in the simulations can be unaffected). As such, the user can adjust the mach
limits for runs to optimize the run time while maintaining sufficient accuracy, possibly
up to mach numbers of 0.1 or more in some cases. Namely, the higher the mach number limits,
the faster the run but the lower the quantitative accuracy. Care should be taken not
to set the mach number limits to be too high such that accuracy suffers to the degree
that renders the simulation results unreliable. One could - for example - do test runs
for a given problem with different mach limits to determine how high the mach number can be
before the results are adversely impacted (ie. when the results depart significantly from
those of lower mach number runs such as at mach <~ 0.01).

The optimal mach limits that minimize runtime while retaining sufficient accuracy will be problem
dependent (eg. plastic rheology dependent etc). As such, each problem requires tests to verify
how high the mach number limits can be set. A good starting point is to set mach number
limits at 0.005 and 0.015 (ie. to maintain the mach number around 0.01), and then to do
runs to test how much higher these limits may be set before the results are significantly
affected. As GTXLBM.py has a high throughput and scaling performance, even runs using
the lower mach limits (around 0.01) will be fast. But if for a large problem when runs
times are large despite the high performance of GTXLBM.py, it may be well worthwhile
to study how high the mach number can be pushed while retaining reliable results, as one
may be able to achieve a 2x to 20x runtime reduction!

--------------------------------------------------------------------------------------------
Domain decomposition
--------------------------------------------------------------------------------------------

GTXLBM.py calculates an optimal domain decomposition for a given run by scanning through
all x, y and z domain sizes (nxi, nyi and nzi) to find the most optimal choice which maximizes
core useage and efficiency of sending and receiving domain edges between cores.

The x axis is the domain that is used for cache blocking in the c code, and x domain scan
limits of around nxi_min ~ 200 and nxi_max ~ 512 has been found to be optimal on the computers
I have run on. However, the optimal x domain size for cache blocking may vary depending on
the chips on the computer being used so it may be worth to do a test on your computer if
you are doing large runs.

The y and z scan limits (nyi_min, nyi_max, nzi_min, nzi_max) can usually be left unspecified
for 2D runs, but for larger runs or 3D runs, the scans may start to take a significant amount
of time. In such cases, users can set the ranges of these scans such that the domain
decomposition can be reduced to seconds or a fraction of a second by a suitable choice of
the y and z scan limits. Below, I explain how this can be done.

One should first assume the x domain size denoted nxi will be in the mid range of the
x scan limits (eg. nxi ~ 300). The y and z scan limits should then be chosen around values
that would give an optimal somain decompisiton, so roughly equal sized y and z domain sizes,
so nyi ~ nzi. Knowing this, one can estimate the optimal y and z domain sizes, and as such,
set the y and z scan limits to be around these values so scan time is dramatically reduced.

To illustrate, consider a 3D problem of size nx = 2501, ny = 2501, nz = 201. As nxi denotes
the domain size in x, there will be around nx/nxi = 2501/300 ~ 8 domains in x, so ndx = 8.
By definition, the product of the number of domains in x, y and z must equal the number of
cores being used, so ndx.ndy.ndz = ncores. So if you wish to run on ncores = 4000 cores,
then there will be ~ nx/300 x domains. Hence, ndx.ndy.ndz = 4000 where ndx=8. Hence, ndy.ndz
= 4000/8 = 500. Since the optimal domain sizes in y and z should be about the same, we have
nyi = ny/ndy = nzi = nz/ndz. We also know that that ny = 2501/201 = 12.44 ~ 12. Hence, we
have ny/ndy = 12.nz/ndy = nz/ndz so ndy = 12 ndz. Hence,
we have ndy.ndz = 500 = 12.ndz^2. Hence, ndz = sqrt(500/12) ~ 6. Hence, ndy ~ 12x6 = 72.
Since the z domain size is given by nzi = nz/ndz, we have nzi ~ 201/6 ~ 33, and similarly,
nyi = ny/nyi = 2501/72 ~ 34. As such, the scan ranges should be set to be from somewhat
below to somewhat above 33. For example, nyi_min = 20, nyi_max = 45, nzi_min = 20,
nzi_max = 45.

--------------------------------------------------------------------------------------------
Choosing the number of cores to run on
--------------------------------------------------------------------------------------------

Reference [5] shows that GTXLBM has comparable or better performance to a highly optimized
open source C++ LBM code (OpenLB). Furthermore, [5] shows linear weak scaling performance
is projected up to 100s of Kcores, with the number depending on the interconnect bandwidth.
Tests to date have demonstrated linear weak scaling to 50 Kcores which is the largest number
of cores I have run on so far. However, strong scaling is another matter.

Namely, if you have a fixed problem size and wish to run on more cores, then the initial
linearly scaling performance with number of cores eventually rolls over as the domain sizes
become smaller and the cost of sending and receiving domain edge data exceeds time of
the calculations within the domains. Usually, this occurs when the domain sizes are of
order tens. As such, while GTXLBM.py will achieve linear scaling and high performance
for large 2D of 3D problems using up to tens to hundreds of K cores, the optimal number
of cores for high throughput and good (linear) speedup may be much less for smaller problems.
For examples, it may be optimal to use 16 cores for small 2D problems, and to use of order
hundreds, thousands, or tens of thousands of cores for a large 2D problems and 3D problems.

============================================================================================
References relating to developing and testing the TLBM for geodynamics research
============================================================================================

[1] Mora, P. & Yuen, D. (2018) Simulation of regimes of convection and plume dynamics by the
    thermal Lattice Boltzmann Method, Physics of the Earth & Planetary Interiors 275, 69-79.
    (describes the thermal LBM and shows its ability to model low to high Ra and Pr)

[2] Mora, P., Morra, G. & Yuen, D. (2020) A concise python implementation of the lattice
    Boltzmann method on HPC for geo-fluid flow, Geophysical J. Int. 220 (1), 682-702.
    (desribes basic parallel python implementation of the LBM, and shows scaling on a cluster. 

[3] Mora, P., Morra, G. & Yuen, D. (2023) Models of plate tectonics with the Lattice Boltzmann
    Method, Artificial Intelligence in the Geosciences, 4, 47–58.
    (describes the TLBM and its implementation for highly nonlinear plastic rheologies,
     and demonstrates that plate-tectonic like behavior emerges from plastic rheology alone, ie.
     without the need for temperature dependence of viscosity. Shows how the physical
     time scale of TLBM runs can be calculated)

[4] Morra, G., Mora, P., Hüttig, C., Tosi, N. & Yuen, D. (2024) A benchmark for finite Prantl
    number convection: comparison of Boltzmann and Navier-Stokes solutions, Geophys. J. Int.,
    accepted with minor revision
    (shows modified finite Pr Blanchanbach benchmarks of the TLBM solutions against two
     finite Pr PDE solvers, and the TLBM at Pr=1000 against the infinite Pr solution)

[5] Mora, P., Morra, G. & Yuen, D. (2025) HPC Thermal Lattice Boltzmann Method implementation
    for exascale computational geodynamics, in preparation for J. Computational Physics.
    (details a highly optimized TLBM CPU cluster code, and its performance & scalability
    to 100's of Kcores)

============================================================================================
References on the Lattice Boltzmann Method
============================================================================================

[6]  Succi, S. (2018) The Lattice Boltzmann Equation: For Complex States of Flowing Matter,
         Oxford University Press, http://dx.doi.org/10.1093/oso/9780199592357.001.0001.

[7] Krüger, T., Kusumaatmaja, H., Kuzmin, A., Shardt, O., Silva, G. & Viggen, E.M. (2017) The Lattice
         Boltzmann Method - Principles and Practice, https://dx.doi.org/10.1007/978-3-319-44649-3.

============================================================================================
