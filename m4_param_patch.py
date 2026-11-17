#!/usr/bin/env python3
"""Add M1-M4 pd/bs to the input-parameter whitelist in get_parameters.py.
Run in the code directory. Idempotent: skips if already applied."""
s = open("get_parameters.py").read()

if 'parameter == "M4_bs"' in s:
    print("already applied")
    raise SystemExit

g = "    global TTG_DI_max, TTG_LP_frac, TTG_MP_frac, TTG_HP_frac, fields_every_n_data\n"
assert s.count(g) == 1, ("global anchor", s.count(g))
s = s.replace(g, g + "    global M1_pd, M2_pd, M3_pd, M4_pd, M1_bs, M2_bs, M3_bs, M4_bs\n")

a = '''            elif parameter == "TTG_DI_max":
                valid = True ; TTG_DI_max = float(value)
'''
assert s.count(a) == 1, ("elif anchor", s.count(a))
ins = a
for p in ["M1_pd", "M2_pd", "M3_pd", "M4_pd", "M1_bs", "M2_bs", "M3_bs", "M4_bs"]:
    ins += '            elif parameter == "%s":\n                valid = True ; %s = float(value)\n' % (p, p)
s = s.replace(a, ins)
open("get_parameters.py", "w").write(s)
print("patched: M1-M4 pd/bs now valid input parameters")
