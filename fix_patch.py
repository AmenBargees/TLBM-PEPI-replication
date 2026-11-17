#!/usr/bin/env python3
"""Hard-M4 extraction cap + downward deposit walk (TLBM_rest_kernel.c)
and physical cont_frac threshold (python). Run from the code directory."""
import shutil, sys

# ---------- C kernel ----------
p = "TLBM_rest_kernel.c"
shutil.copy(p, p + ".bak_m4fix")
s = open(p).read()

a1 = "                        M_ext_bs[y][z][x] = M_ext_bs_LP_TTG[y][z][x] + M_ext_bs_MP_TTG[y][z][x] + M_ext_bs_HP_TTG[y][z][x];\n"
assert s.count(a1) == 1, ("a1", s.count(a1))
ins1 = a1 + """
                        /* HARD cumulative cap (Gunawardana M4 semantics): once this
                           cell has delivered M4 of TTG it is infertile - nothing
                           more is extracted or deposited from here. */
                        tmp = M4_bs - cumulative_extracted_bs[y][z][x];
                        if (tmp < 0) {tmp = 0;}
                        if (M_ext_bs[y][z][x] > tmp) {
                            FLOAT _s = (M_ext_bs[y][z][x] > 0) ? tmp / M_ext_bs[y][z][x] : 0;
                            M_ext_bs_LP_TTG[y][z][x] *= _s;
                            M_ext_bs_MP_TTG[y][z][x] *= _s;
                            M_ext_bs_HP_TTG[y][z][x] *= _s;
                            M_ext_bs[y][z][x]         = tmp;
                        }
"""
s = s.replace(a1, ins1)

a2 = "                        if (M_ext_pd[y][z][x] > tmp) {M_ext_pd[y][z][x] = tmp;}\n"
assert s.count(a2) == 1, ("a2", s.count(a2))
ins2 = a2 + """                        /* HARD cumulative cap for peridotite (same semantics). */
                        tmp = M4_pd - cumulative_extracted_pd[y][z][x];
                        if (tmp < 0) {tmp = 0;}
                        if (M_ext_pd[y][z][x] > tmp) {M_ext_pd[y][z][x] = tmp;}
"""
s = s.replace(a2, ins2)

old_dep = """                        if (strcmp(melt_position,"surface") == 0 && boundary_u) {
                            int   zd_dep = zu + 1;
                            FLOAT take   = M_ext_bs[y][z][x];
                            if (take > rho_b[y][zd_dep][x]) take = rho_b[y][zd_dep][x];
                            if (take < 0)                    take = 0;
                            rho_b[y][zd_dep][x] -= take;
                            rho_r[y][zd_dep][x] += take;
                            for (int _aa = 0; _aa < na; _aa++) {
                                f_b[y][zd_dep][_aa][x] -= take * w[_aa];
                                f_r[y][zd_dep][_aa][x] += take * w[_aa];
                            }
                        }
"""
assert s.count(old_dep) == 1, ("dep", s.count(old_dep))
new_dep = """                        if (strcmp(melt_position,"surface") == 0 && boundary_u) {
                            /* Deposit walks DOWN from zu+1 until it finds basalt to
                               convert (crustal thickening) instead of silently
                               dropping melt on a saturated surface. Keeps the red
                               budget equal to the bookkept extraction. */
                            FLOAT rem = M_ext_bs[y][z][x];
                            int   zd_dep;
                            for (zd_dep = zu + 1; zd_dep <= zu + 40 && rem > 0; zd_dep++) {
                                FLOAT take = rem;
                                if (take > rho_b[y][zd_dep][x]) take = rho_b[y][zd_dep][x];
                                if (take < 0)                    take = 0;
                                rho_b[y][zd_dep][x] -= take;
                                rho_r[y][zd_dep][x] += take;
                                for (int _aa = 0; _aa < na; _aa++) {
                                    f_b[y][zd_dep][_aa][x] -= take * w[_aa];
                                    f_r[y][zd_dep][_aa][x] += take * w[_aa];
                                }
                                rem -= take;
                            }
                        }
"""
s = s.replace(old_dep, new_dep)
open(p, "w").write(s)
print("C kernel patched")

# ---------- python cont_threshold ----------
for pf in sys.argv[1:]:
    t = open(pf).read()
    old = "cont_threshold     = 0\n"
    assert t.count(old) == 1, (pf, t.count(old))
    t = t.replace(old, "cont_threshold     = 0.5   # real crust, not a film\n")
    open(pf, "w").write(t)
    print("threshold patched:", pf)
