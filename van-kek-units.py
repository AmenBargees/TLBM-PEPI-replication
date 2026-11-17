import numpy as np

g_p = 9.8
g_l = 1.952e-5
nu_l = 0.7889
ka_l = 0.0007889
Pr_l  = 1e3
rho_p = 3200
eta_p = 1e21
nu_p  = eta_p/rho_p
dt_p_file = 7.889e8
dx_p  = 1e3
nx    = 501
L_p   = 0.9142
L_l   = int((nx-1)*L_p+0.5)

Lx    = (nx-1)*dx_p
dt_p  = g_l*nu_p/(g_p*nu_l*dx_p) 
t_p   = 5*1e6*365*24*3600 ## 5 My similar to t=20
nt    = t_p/dt_p

print("-----------------")
print("rk2-van-kek")
print("-----------------")
print("dt_p = %8.3g" %(dt_p))
print("L_p = %8.3g" %(L_p))
print("L_l = %8.3g" %(L_l))
print("Pr_l = %8.3g" %(Pr_l))
print("ka_l = %8.3g" %(ka_l))
print("-----------------")
print("nu_l = %8.3g" %(nu_l))
print("g_l  = %8.3g" %(g_l))
print("-----------------")
print("-----------------")


print("dt_p = %8.3g ?= %8.3g = from out file" %(dt_p,dt_p_file))
nt    = 200000
print("nt from file = ",nt)
dt_p = t_p/(nt-1)
print("dt_p = t_p/(nt-1) = %8.3g " %(dt_p))

Lx_kek  = 1
t_kek   = 2000
eta_kek = 100
rho_kek = 1000
nu_kek  = eta_kek/rho_kek

factor_t  = t_p /t_kek
factor_x  = Lx  /Lx_kek
factor_nu = nu_p/nu_kek

print("factor_t  = %8.3g" %(factor_t))
print("factor_x  = %8.3g" %(factor_x))
print("factor_nu = %8.3g" %(factor_nu))

dt_kek = t_kek/(nt-1)
print("dt_kek = %8.3g" %(dt_kek))
print("nt            = ",nt)
print("(nt-1)*dt_kek = %8.3g" %((nt-1)*dt_kek))


eta_p = eta_kek
rho_p = rho_kek
dt_p  = dt_kek*0.3
nu_p  = eta_p/rho_p
ka_p  = 1e-6
C_p   = rho_p*3e-5
C_l   = C_p
Ra_p  = g_p*C_p*L_p**3/(nu_p*ka_p)
Pr_p  = nu_p/ka_p
Ra_l  = Ra_p
ka_l  = dt_p*L_l**2*ka_p/L_p**2
nu_l  = ka_l*Pr_l
g_l   = nu_l**2*Ra_l/(C_l*L_l**3*Pr_l)

print("-----------------")
print("dt_p = %8.3g" %(dt_p))
print("L_p = %8.3g" %(L_p))
print("L_l = %8.3g" %(L_l))
print("Pr_l = %8.3g" %(Pr_l))
print("ka_l = %8.3g" %(ka_l))
print("-----------------")
print("nu_l = %8.3g" %(nu_l))
print("g_l  = %8.3g" %(g_l))
print("-----------------")

