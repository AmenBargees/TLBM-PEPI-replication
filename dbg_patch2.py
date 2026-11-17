import re, shutil, sys
import sys; p = sys.argv[1]
shutil.copy(p, p + ".bak_dbg")
s = open(p).read()

helper = '''
def _dbg_redsum(tag):
    _R = get_field(rho_r)
    if rank == 0:
        print("DBGRED %-6s it=%d sum_r=%.9e" % (tag, it, float(np.sum(_R, dtype=np.float64))), flush=True)

'''
anchor = "while int(t-dt) <= nt:"
assert s.count(anchor) == 1
s = s.replace(anchor, helper + anchor)

# gate + preK after the it = int(...) line
m = re.search(r"\n(    it = int\(t\*dt_dt0\)[^\n]*)\n", s)
assert m, "it anchor not found"
gate = ("\n" + m.group(1) + "\n"
        "    _dbg_on = comm.bcast(bool((2000 <= it <= 2020) or (2245 <= it <= 2262) or (30000 <= it <= 30030) or (60000 <= it <= 60030) or (100000 <= it <= 100030) or ((it % 2000) < 3)), root=0)\n"
        "    if _dbg_on: _dbg_redsum(\"preK\")\n")
s = s[:m.start()] + gate + s[m.end():]

# postK after the RKTLBM_rest call's swaps (first totalTimeCalcs after gatev close)
k = s.find("np_float(lp_gatev), np_float(mp_gatev), np_float(hp_gatev))")
assert k > 0
j = s.find("totalTimeCalcs += (time.time()-t0)", k)
assert j > 0
eol = s.find("\n", j)
s = s[:eol] + "\n                if _dbg_on: _dbg_redsum(\"postK\")" + s[eol:]

open(p, "w").write(s)
print("patched OK")
