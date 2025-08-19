import sys
import cmath
import math
import os
import h5py
import matplotlib.pyplot as plt   # plots
import numpy as np
from scipy.interpolate import griddata

from liblibra_core import *
import util.libutil as comn
from libra_py import units
from libra_py import data_conv
import libra_py.dynamics.tsh.compute as tsh_dynamics
import libra_py.data_savers as data_savers

from plot_utils import set_style, add_hbar, add_abar

# Add the parent directory to sys.path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir)

from kcrpmd_utils.kcrpmdtst import KcrpmdTst
from kcrpmd_utils.kcrpmdmodel import gen_kcrpmd_bath_params, get_ABC, kcrpmd_system_bath

set_style()

calc_dirs = [d for d in os.listdir('../') if d.startswith('_sys_')]
_sys_3_method_1_K0hwp1 = sorted([k for k in calc_dirs if k.startswith("_sys_3_method_1") and "_hw_1" in k and "_K0_1.00e-10" in k], key=lambda s: float(s.split('_')[10]))
_sys_3_method_1_K0hwn1 = sorted([k for k in calc_dirs if k.startswith("_sys_3_method_1") and "_hw_-1" in k and "_K0_1.00e-10" in k], key=lambda s: float(s.split('_')[10]))
_sys_3_method_1_Khwp1 = sorted([k for k in calc_dirs if k.startswith("_sys_3_method_1") and "_hw_1" in k and not "_K0_1.00e-10" in k], key=lambda s: float(s.split('_')[10]))
_sys_3_method_1_Khwn1 = sorted([k for k in calc_dirs if k.startswith("_sys_3_method_1") and "_hw_-1" in k and not "_K0_1.00e-10" in k], key=lambda s: float(s.split('_')[10]))
_sys_3_method_3_Khwp1 = sorted([k for k in calc_dirs if k.startswith("_sys_3_method_3") and "_hw_1" in k], key=lambda s: float(s.split('_')[8]))
_sys_3_method_3_Khwn1 = sorted([k for k in calc_dirs if k.startswith("_sys_3_method_3") and "_hw_-1" in k], key=lambda s: float(s.split('_')[8]))

with open("../" + _sys_3_method_1_K0hwp1[0] + "/_model_params.txt") as f:
    model_params = eval(f.read())

with open("../" + _sys_3_method_1_K0hwp1[0] + "/_control_params_dynamics.txt") as f:
    control_params = eval(f.read())

# ======= Pull in all the rate data =======
beta = units.hartree / (units.boltzmann * control_params["Temperature"])

leps_arr = np.array([key.split('_')[10] for key in _sys_3_method_1_K0hwp1], dtype=float)[:]
kBO0_arr = np.zeros(leps_arr.shape)
kBO0_se_arr = np.zeros(leps_arr.shape)

kGR_arr = np.zeros(leps_arr.shape)
kBO_arr = np.zeros(leps_arr.shape)
kBO_se_arr = np.zeros(leps_arr.shape)
kIF_arr = np.zeros(leps_arr.shape)
kIF_se_arr = np.zeros(leps_arr.shape)
knew_arr = np.zeros(leps_arr.shape)
knew_se_arr = np.zeros(leps_arr.shape)

for i, d in enumerate(_sys_3_method_1_K0hwp1):
    #kGR0_arr[i] = np.loadtxt("../" + _sys_3_method_1_K0hwp1[i] + "/tst_data/kGR_full.txt")
    Phw0_p1 = np.loadtxt("../" + _sys_3_method_1_K0hwp1[i] + "/tst_data/Phw.txt")
    ktst0_p1 = np.loadtxt("../" + _sys_3_method_1_K0hwp1[i] + "/tst_data/ktsts.txt")
    #kappa0_p1 = np.loadtxt("../" + _sys_3_method_1_K0hwp1[i] + "/kappa_data/kappa_avg.txt")[-1]
    kappa0_p1 = 1.0 
    #kappa0_se_p1 = np.loadtxt("../" + _sys_3_method_1_K0hwp1[i] + "/kappa_data/kappa_se.txt")[-1]
    kappa0_se_p1 = 0.0
    Phw0_n1 = np.loadtxt("../" + _sys_3_method_1_K0hwn1[i] + "/tst_data/Phw.txt")
    ktst0_n1 = np.loadtxt("../" + _sys_3_method_1_K0hwn1[i] + "/tst_data/ktsts.txt")
    #kappa0_n1 = np.loadtxt("../" + _sys_3_method_1_K0hwn1[i] + "/kappa_data/kappa_avg.txt")[-1]
    kappa0_n1 = 1.0 
    #kappa0_se_n1 = np.loadtxt("../" + _sys_3_method_1_K0hwn1[i] + "/kappa_data/kappa_se.txt")[-1]
    kappa0_se_n1 = 0.0
    kBO0_p1 = Phw0_p1 * ktst0_p1 * kappa0_p1
    kBO0_n1 = Phw0_n1 * ktst0_n1 * kappa0_n1
    kBO0_se_p1 = Phw0_p1 * ktst0_p1 * kappa0_se_p1
    kBO0_se_n1 = Phw0_n1 * ktst0_n1 * kappa0_se_n1
    kBO0_arr[i] = kBO0_p1 + kBO0_n1
    kBO0_se_arr[i] = np.sqrt(kBO0_se_p1**2 + kBO0_se_n1**2)
    #print(kBO0_p1, kBO0_se_p1, kBO0_n1, kBO0_se_n1)
    #######
    #kBO0_arr[i] = Phw0_p1 * ktst0_p1 + Phw0_n1 * ktst0_n1
    #kBO0_se_arr[i] = 0.0 
    #######

#print(_sys_3_method_1_K0hwp1)
#print(_sys_3_method_1_Khwp1)
#print(_sys_3_method_1_K0hwn1)
#print(_sys_3_method_1_Khwn1)
#exit()

for i, d in enumerate(_sys_3_method_1_Khwp1):
    kGR_arr[i] = np.loadtxt("../" + _sys_3_method_1_Khwp1[i] + "/tst_data/kGR_full.txt")
    Phw_p1 = np.loadtxt("../" + _sys_3_method_1_Khwp1[i] + "/tst_data/Phw.txt")
    ktst_p1 = np.loadtxt("../" + _sys_3_method_1_Khwp1[i] + "/tst_data/ktsts.txt")
    kappa_p1 = np.loadtxt("../" + _sys_3_method_1_Khwp1[i] + "/kappa_data/kappa_avg.txt")[-1]
    kappa_se_p1 = np.loadtxt("../" + _sys_3_method_1_Khwp1[i] + "/kappa_data/kappa_se.txt")[-1]
    Phw_n1 = np.loadtxt("../" + _sys_3_method_1_Khwn1[i] + "/tst_data/Phw.txt")
    ktst_n1 = np.loadtxt("../" + _sys_3_method_1_Khwn1[i] + "/tst_data/ktsts.txt")
    kappa_n1 = np.loadtxt("../" + _sys_3_method_1_Khwn1[i] + "/kappa_data/kappa_avg.txt")[-1]
    kappa_se_n1 = np.loadtxt("../" + _sys_3_method_1_Khwn1[i] + "/kappa_data/kappa_se.txt")[-1]
    kBO_p1 = Phw_p1 * ktst_p1 * kappa_p1
    kBO_n1 = Phw_n1 * ktst_n1 * kappa_n1
    kBO_se_p1 = Phw_p1 * ktst_p1 * kappa_se_p1
    kBO_se_n1 = Phw_n1 * ktst_n1 * kappa_se_n1
    kBO_arr[i] = kBO0_p1 + kBO_n1
    kBO_se_arr[i] = np.sqrt(kBO_se_p1**2 + kBO_se_n1**2)

for i, d in enumerate(_sys_3_method_3_Khwp1):
    print(d)
    Phw_p1 = np.loadtxt("../" + _sys_3_method_3_Khwp1[i] + "/tst_data/Phw.txt")
    ktst_p1 = np.loadtxt("../" + _sys_3_method_3_Khwp1[i] + "/tst_data/ktsts.txt")
    kappa_p1 = np.loadtxt("../" + _sys_3_method_3_Khwp1[i] + "/kappa_data/kappa_avg.txt")[-1]
    kappa_se_p1 = np.loadtxt("../" + _sys_3_method_3_Khwp1[i] + "/kappa_data/kappa_se.txt")[-1]

    Phw_n1 = np.loadtxt("../" + _sys_3_method_3_Khwn1[i] + "/tst_data/Phw.txt")
    ktst_n1 = np.loadtxt("../" + _sys_3_method_3_Khwn1[i] + "/tst_data/ktsty.txt")
    kappa_n1 = np.loadtxt("../" + _sys_3_method_3_Khwn1[i] + "/kappa_data/kappa_avg.txt")[-1]
    kappa_se_n1 = np.loadtxt("../" + _sys_3_method_3_Khwn1[i] + "/kappa_data/kappa_se.txt")[-1]

    knew_p1 = Phw_p1 * ktst_p1 * kappa_p1
    knew_n1 = Phw_n1 * ktst_n1 * kappa_n1
    knew_se_p1 = Phw_p1 * ktst_p1 * kappa_se_p1
    knew_se_n1 = Phw_n1 * ktst_n1 * kappa_se_n1
    knew_arr[i] = knew_p1 + knew_n1
    knew_se_arr[i] = np.sqrt(knew_se_p1**2 + knew_se_n1**2)

print(kGR_arr)
print(kBO0_arr)
print(kBO0_se_arr)
print(kBO_arr)
print(kBO_se_arr)
print(knew_arr)
print(knew_se_arr)
exit()



#for i, d in enumerate(_sys_3_method_1_K0hwp1):

    #ktst0_p1 = np.loadtxt("../" + _sys_3_method_1_K0hwp1[0] + "/tst_data/ktsts.txt")
    #kappa0_avg_p1 = np.loadtxt("../" + _sys_3_method_1_K0hwp1[0] + "/kappa_data/kappa_avg.txt")[-1]
    #kappa0_se_p1 = np.loadtxt("../" + _sys_3_method_1_K0hwp1[0] + "/kappa_data/kappa_se.txt")[-1]
    #kBO0_p1 = ktst0_p1 * kappa0_avg_p1
    #kBO0_se_p1 = ktst0_p1 * kappa0_se_p1
    #Phw_n1 = np.loadtxt("../" + _sys_3_method_1_K0hwn1[0] + "/tst_data/Phw.txt")
    #ktst0_n1 = np.loadtxt("../" + _sys_3_method_1_K0hwn1[0] + "/tst_data/ktsts.txt")
    #kappa0_avg_n1 = np.loadtxt("../" + _sys_3_method_1_K0hwn1[0] + "/kappa_data/kappa_avg.txt")[-1]
    #kappa0_se_n1 = np.loadtxt("../" + _sys_3_method_1_K0hwn1[0] + "/kappa_data/kappa_se.txt")[-1]
    #kBO0_n1 = ktst0_n1 * kappa0_avg_n1
    #kBO0_se_n1 = ktst0_n1 * kappa0_se_n1


#K0_arr = np.array([key.split('_')[8] for key in _sys_3_method_1], dtype=float)[1:]


for i, d in enumerate(_sys_3_method_1[1:]):
    kGR_arr[i] = np.loadtxt("../" + _sys_3_method_1[i+1] + "/tst_data/kGR.txt")
    kBO_arr[i] = np.loadtxt("../" + _sys_3_method_1[i+1] + "/tst_data/ktsts.txt")
    kBO_se_arr[i] = kBO_arr[i] * np.loadtxt("../" + _sys_3_method_1[i+1] + "/kappa_data/kappa_se.txt")[-1]
    kBO_arr[i] *= np.loadtxt("../" + _sys_3_method_1[i+1] + "/kappa_data/kappa_avg.txt")[-1]
    kIF_arr[i] = kGR_arr[i] * kBO_arr[i] / (kGR_arr[i] + kBO0)

for i in range(4):
    knew_arr[i] = np.loadtxt("../" + _sys_3_method_3_fix_y[i] + "/tst_data/ktsty.txt")
    knew_se_arr[i] = knew_arr[i] * np.loadtxt("../" + _sys_3_method_3_fix_y[i] + "/kappa_data/kappa_se.txt")[-1]
    knew_arr[i] *= np.loadtxt("../" + _sys_3_method_3_fix_y[i] + "/kappa_data/kappa_avg.txt")[-1]
for i in range(4,10):
    knew_arr[i] = np.loadtxt("../" + _sys_3_method_3_fix_s[i] + "/tst_data/ktsts.txt")
    knew_se_arr[i] = knew_arr[i] * np.loadtxt("../" + _sys_3_method_3_fix_s[i] + "/kappa_data/kappa_se.txt")[-1]
    knew_arr[i] *= np.loadtxt("../" + _sys_3_method_3_fix_s[i] + "/kappa_data/kappa_avg.txt")[-1]

fig, ax = plt.subplots()
ax.errorbar(np.log10(beta * K0_arr), np.log10(kGR_arr), kGR_arr*0, fmt='o-', markersize=3, linewidth=1, color='r', label=r'$k_\mathrm{GR}$')
ax.errorbar(np.log10(beta * K0_arr), np.log10(kBO_arr), kBO_se_arr / (kBO_arr * np.log(10)), fmt='o-', markersize=3, linewidth=1, color='b', label=r'$k_\mathrm{BO}$')
ax.errorbar(np.log10(beta * K0_arr), np.log10(kIF_arr), kIF_se_arr, fmt='^-', markersize=3, linewidth=1, color='k', label=r'$k_\mathrm{IF}$')
ax.errorbar(np.log10(beta * K0_arr), np.log10(knew_arr), knew_se_arr / (knew_arr * np.log(10)), fmt='*-', markersize=3, linewidth=1, color='gold', label=r'$k_\mathrm{new}$')
#ax.set_xlim(-0.3,1.3)
#ax.set_ylim(-8.0,5.0)
#ax.set_xticks([0.0, 1.0])
#ax.set_yticks([0])
#ax.set_xticklabels([r'$0$', r'$q_0$'])
ax.set_xlabel(r"log(β$K_0$)", fontsize = 15)
ax.set_ylabel(r"log($k_{\mathrm{ET}}$)", fontsize = 15)
ax.set_title("")
ax.legend(loc='upper left')

plt.tight_layout()
#plt.subplots_adjust(left=0.1, right=0.98, top=0.98, bottom=0.18)
plt.savefig('fig5.png')



