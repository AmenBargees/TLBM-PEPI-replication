
import numpy as np
import matplotlib
import matplotlib.pyplot as plt

matplotlib.rcParams['text.usetex'] = False

cores_huge       = np.array([ 75032, 99830, 148808, 199980, 299710],dtype=int)
mflups_core_huge = np.array([  6.24,  6.08,   6.00,   4.93,   4.51])
# Runs/nodes                   391a,   521,   782a,   1042a,  1563,   2084

mflups_huge = cores_huge * mflups_core_huge
Kcores = cores_huge/1000

nx = 13525
ny = 13525
nz = 967
n = float(nx)*float(ny)*float(nz)
print("n_huge = %6.3g" %(n))

print("n_575  = %6.3g" %(float(575))**3)
print("n_1449 = %6.3g" %(float(1449))**3)

Tflups_huge = mflups_huge/1.e6

Tlattice_sites = n/1.e12

secs_per_timestep = Tlattice_sites/Tflups_huge

print("Tflups_huge            = ",Tflups_huge)
print("secs_per_timestep_huge = ",secs_per_timestep)

mach = 0.01
nx_wilson = 2*np.pi*(nz-1)/2
nt_wilson = nx_wilson/mach
days_wilson = nt_wilson*secs_per_timestep/3600/24
print("nt_wilson /1000       = ",int(nt_wilson/1.e3+1))
print("days per wilson cycle = ",days_wilson)

print("Shaheen III")
print("cores_huge = ",cores_huge)
print("nodes_huge = ",np.int32(cores_huge/192+0.999))
print("mflups_huge = ",mflups_huge)

# OpenLB Ice Lake

nodes_OLB = np.array(           [   1,   8,  16,  24,  32,  40,  48,  56,  64,  72,  80,  88,  96, 104, 112, 120, 128])
gflups_node_OLB_575 =  np.array([1.25,0.82,0.74,0.60,0.59,0.56,0.54,0.52,0.50,0.48,0.47,0.46,0.39,0.38,0.40,0.40,0.40])
gflups_node_OLB_1449 = np.array([ 1.3,1.15,0.92,0.92,0.90,0.88,0.87,0.86,0.93,0.85,0.84,0.83,0.82,0.85,0.83,0.82,0.82])
cores_OLB = nodes_OLB * 76
mflups_core_OLB_575 = gflups_node_OLB_575*1000/(2*76)
mflups_core_OLB_1449 = gflups_node_OLB_1449*1000/(2*76)
#plt.plot(cores_OLB,mflups_core_OLB_575,label="575^3")
#plt.plot(cores_OLB,mflups_core_OLB_1449,label="1449^3")
#plt.legend()
#plt.plot(nodes_OLB,gflups_node_OLB_575,label="575^3")
#plt.plot(nodes_OLB,gflups_node_OLB_1449,label="1449^3")
#plt.legend()
#plt.show()

# Shaheen III OpenLB

cores_575       = np.array([ 192*1, 192*32, 192*64, 192*128, 192*256])
mflups_core_575 = np.array([   4.8,   3.93,   3.75,    3.37,    2.7])
cores_575       = np.array([ 192*1, 192*2, 192*4, 192*8, 192*16, 192*32, 192*64, 192*128, 192*256])
mflups_core_575 = np.array([   4.8,  4.52,  4.29,  4.13,   3.88,   3.77,   3.67,    3.37,    2.7])
mflups_575      = mflups_core_575 * cores_575
Kcores_575      = cores_575/1000

#cores_913       = np.array([ 192*1, 192*2, 192*4, 192*8, 192*16, 192*32, 192*64, 192*128, 192*256])
#mflups_core_913 = np.array([  4.7,  4.66,  4.63,  4.54,   4.2,   4.14,   3.97  ,    3.61,    3.58])
cores_913       = np.array([ 192*1, 192*32, 192*64, 192*128, 192*256])
mflups_core_913 = np.array([  4.7,  4.14,   3.97  ,    3.61,    3.58])
mflups_913      = mflups_core_913 * cores_913
Kcores_913      = cores_913/1000

cores_1449       = np.array([ 192*4, 192*32, 192*64, 192*128, 192*256])
mflups_core_1449 = np.array([ 4.59,  4.39,   4.39  ,    4.00,    3.91])
cores_1449       = np.array([ 192*4, 192*8, 192*16, 192*32, 192*64, 192*128, 192*256])
mflups_core_1449 = np.array([ 4.59,  4.51,   4.43,   4.39,   4.39  ,    4.00,    3.91])
mflups_1449      = mflups_core_1449 * cores_1449
Kcores_1449      = cores_1449/1000

cores_per_node   = 2*38

nodes_OpenLB  = np.int32(cores_575/cores_per_node+0.5)
nodes_Shaheen = np.int32(cores_575/192+0.999)

gflups_node_575  = 2*mflups_core_575*cores_per_node/1000
gflups_node_913  = 2*mflups_core_913*cores_per_node/1000
gflups_node_1449 = 2*mflups_core_1449*cores_per_node/1000

print("cores = ",cores_575)
print("nodes_OpenLB  = ",nodes_OpenLB)
print("gflups_node   = ",gflups_node_575)

fig, ax = plt.subplots(figsize=(7,5))
ax.set_xscale("log")
axR = ax.twinx()

axT = ax.twiny()
axT.set_xscale("log")

ax.grid(linestyle="-",which="both")
ax.grid(which='minor', alpha=0.4)
ax.grid(which='major', alpha=1)

nodes_OLB_benchmark = np.array([1,128],dtype=float)
mflups_OLB_benchmark = np.array([0,7],dtype=float)

axT.fill_between(nodes_OLB_benchmark, 7, where=mflups_OLB_benchmark>=0,facecolor='grey', alpha=0.1)

ax.plot(cores_575,mflups_core_575,"-o",color="r"     ,label="575^3  (TLBM)")
ax.plot(cores_OLB,mflups_core_OLB_575,"--x",color="r" ,label="575^3  (OpenLB)")
ax.plot(cores_1449,mflups_core_1449,"-o",color="g"   ,label="1449^3 (TLBM)")
ax.plot(cores_OLB,mflups_core_OLB_1449,"--x",color="g",label="1449^3 (OpenLB)")
ax.plot(cores_huge,mflups_core_huge,"-o",color="b"   ,label="13,525^2 x 967 (TLBM)")
#ax.plot(cores_913,mflups_core_913,"-o",color="b",label="913")


#ax.axvline(192,linestyle="--",label="192 cores",color="r")
#ax.axvline(128*38*2,linestyle="--",label="  9,782 cores = 128 HoreKa eqiv. nodes",color="g")
#ax.axvline(256*192,linestyle="--" ,label="49,192 cores",color="b")

#axR.plot(cores_1449,gflups_node_1449,"-o",color="r",label="1449 (HoreKa)")
#axR.plot(cores_913,gflups_node_913,"-o",color="r",label="913 (HoreKa)")
#axR.plot(cores_575,gflups_node_575,"--o",color="r",label="575 (HoreKa)")
xmax = 1.e6
ax.set_xlim(100,xmax)
ax.set_ylim(0,7)
axT.set_xlim(100/(2*38),xmax/(2*38))
axR.set_ylim(0,2*7*cores_per_node/1000)
ax.set_xlabel("Cores") # ,color="blue")
ax.set_ylabel("TLBM Mflups/core") # ,color="blue")
ax.tick_params(axis='x') # ,labelcolor="blue")
ax.tick_params(axis='y') # ,labelcolor="blue")
axT.set_xlabel('HoreKa equiv. nodes ') # ,color="red")
axR.set_ylabel('(LBM Gflups)/(HoreKa eqiv. node)') # ,color="red")
axT.tick_params(axis='x') # ,labelcolor="red")
axR.tick_params(axis='y') # ,labelcolor="red")
ax.legend(loc="lower right",framealpha=1,facecolor="white")
#ax.legend(loc="best",framealpha=1)
plt.savefig("OpenLB_compare2.png")
plt.close()
exit()
#plt.show()

plt.grid()
plt.plot(Kcores_575,mflups_575/1.e3,"-o",color="r",label="575")
plt.plot(Kcores_913,mflups_913/1.e3,"-o",color="g",label="913")
plt.plot(Kcores_1449,mflups_1449/1.e3,"-o",color="b",label="1449")
#plt.xscale("log")
#plt.yscale("log")
plt.xlabel("K cores")
plt.ylabel("Gflups")
plt.xlim(0,25.0)
plt.ylim(0,100)
plt.legend(loc="upper left")
plt.show()

