#!/usr/bin/env python3
"""Canonical vK97 RT initial condition for vk_b5L on the GTXLBM_BENCH stack:
light (red) layer of thickness 0.2 at the bottom, interface at
y = 0.2 + 0.02*cos(pi*x/Lx). z row 0 = top. rho1 = rho1_p/rho2_p (auto-derived
by the code), rho2 = 1."""
import numpy as np, os
nx, nz = 915, 1001
Lx, Ly = 0.9142, 1.0
rho1 = 1.0 ; rho2 = 1.0   # V26 convention: identity fields, uniform total density
x = (np.arange(nx)+0.5)/nx*Lx
y_int = 0.2 + 0.02*np.cos(np.pi*x/Lx)          # height of interface above base
z = np.arange(nz)                               # 0 = top
h_above_base = (nz-1-z)/(nz-1)*Ly               # height of each row above base
red = (h_above_base[:,None] < y_int[None,:])    # True in light bottom layer
RHO_r = np.where(red, rho1, 0.0)[None,:,:].astype(np.float32)
RHO_b = np.where(red, 0.0, rho2)[None,:,:].astype(np.float32)
out = os.path.expanduser("~/src/GTXLBM_GAUGE/projects/RT/vk_b5L/")
os.makedirs(out, exist_ok=True)
np.save(out+"vk_rho_r.npy", RHO_r)
np.save(out+"vk_rho_b.npy", RHO_b)
print("saved IC:", RHO_r.shape, "red fraction %.4f" % (red.mean()))
