#!/usr/bin/env python3
"""ICs for van-keken-miscible-2x on the GTXLBM_BENCH stack (archean2 route):
immiscible vK RT in the colors (red light 1000, blue 2000, canonical 0.2 layer,
0.02 cos interface) + passive miscible circular inclusion in L (Peter's V26
van-keken-miscible blob: center x0=3(nx-1)/4, z0=(nz-1)/4, radius nx/10)."""
import numpy as np, os
nx, nz, ny = 915, 1001, 1
Lx = 0.9142
rho1 = 1.0 ; rho2 = 1.0   # V26 convention: identity fields, uniform total density           # lattice densities (auto-derive convention)
x = (np.arange(nx)+0.5)/nx*Lx
y_int = 0.2 + 0.02*np.cos(np.pi*x/Lx)
z = np.arange(nz)
h = (nz-1-z)/(nz-1)                          # height above base (z=0 is top row)
red = (h[:,None] < y_int[None,:])
RHO_r = np.where(red, rho1, 0.0)[None].astype(np.float32)
RHO_b = np.where(red, 0.0, rho2)[None].astype(np.float32)
# miscible inclusion in L (exact V26 geometry, same index orientation)
x0, z0, r0 = 3*(nx-1)//4, (nz-1)//4, nx//10
X, Z = np.meshgrid(np.arange(nx), np.arange(nz))
L = (np.sqrt((X-x0)**2 + (Z-z0)**2) < r0).astype(np.float32)[None]
out = os.path.expanduser("~/src/GTXLBM_GAUGE/projects/RT/van-keken-miscible-2x/")
os.makedirs(out, exist_ok=True)
np.save(out+"vkm_rho_r.npy", RHO_r)
np.save(out+"vkm_rho_b.npy", RHO_b)
np.save(out+"vkm_L.npy", L)
print("saved: red frac %.4f, L blob frac %.4f" % (red.mean(), L.mean()))
