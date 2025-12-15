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

# Add the parent directory to sys.path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir)

from kcrpmd_utils.kcrpmdtst import KcrpmdTst

calc_dirs = [d for d in os.listdir('../') if d.startswith('_sys_')]
_sys_1_method_1 = sorted([k for k in calc_dirs if k.startswith("_sys_1_method_1")], key=lambda s: float(s.split('_')[8]))
_sys_1_method_2_fix_y = sorted([k for k in calc_dirs if k.startswith("_sys_1_method_2_fix_y")], key=lambda s: float(s.split('_')[8]))
_sys_1_method_2_fix_s = sorted([k for k in calc_dirs if k.startswith("_sys_1_method_2_fix_s")], key=lambda s: float(s.split('_')[8]))
_sys_1_method_3_fix_y = sorted([k for k in calc_dirs if k.startswith("_sys_1_method_3_fix_y")], key=lambda s: float(s.split('_')[8]))
_sys_1_method_3_fix_s = sorted([k for k in calc_dirs if k.startswith("_sys_1_method_3_fix_s")], key=lambda s: float(s.split('_')[8]))

with open("../" + _sys_1_method_1[0] + "/_model_params.txt") as f:
    model_params = eval(f.read())

with open("../" + _sys_1_method_1[0] + "/_control_params.txt") as f:
    control_params = eval(f.read())

# ======= Pull in all the rate data =======
beta = units.hartree / (units.boltzmann * control_params["Temperature"])

ktst0 = np.loadtxt("../" + _sys_1_method_1[0] + "/tst_data/ktsts.txt")
kappa0_avg = np.loadtxt("../" + _sys_1_method_1[0] + "/kappa_data/kappa_avg.txt")[-1]
kappa0_se = np.loadtxt("../" + _sys_1_method_1[0] + "/kappa_data/kappa_se.txt")[-1]
kBO0 = ktst0 * kappa0_avg
kBO0_se = ktst0 * kappa0_se

K0_array = np.array([key.split('_')[8] for key in _sys_1_method_1], dtype=float)[1:]
kGR_array = np.zeros(K0_array.shape)
kBO_array = np.zeros(K0_array.shape)
kBO_se_array = np.zeros(K0_array.shape)
kIF_array = np.zeros(K0_array.shape)
kIF_se_array = np.zeros(K0_array.shape)
kold_array = np.zeros(K0_array.shape)
kold_se_array = np.zeros(K0_array.shape)
knew_array = np.zeros(K0_array.shape)
knew_se_array = np.zeros(K0_array.shape)

#print(_sys_1_method_1)
#print(len(_sys_1_method_2))
#print(_sys_1_method_2_fix_y)
#print(_sys_1_method_3)

for i, d in enumerate(_sys_1_method_1[1:]):
    kGR_array[i] = np.loadtxt("../" + _sys_1_method_1[i+1] + "/tst_data/kGR.txt")
    kBO_array[i] = np.loadtxt("../" + _sys_1_method_1[i+1] + "/tst_data/ktsts.txt")
    kBO_array[i] *= np.loadtxt("../" + _sys_1_method_1[i+1] + "/kappa_data/kappa_avg.txt")[-1]
    kBO_se_array[i] = kBO_array[i] * np.loadtxt("../" + _sys_1_method_1[i+1] + "/kappa_data/kappa_se.txt")[-1]
    kIF_array[i] = kGR_array[i] * kBO_array[i] / (kGR_array[i] + kBO0)
    kIF_se_array[i] = kBO_se_array[i]

for i in range(5):
    kold_array[i] = np.loadtxt("../" + _sys_1_method_2_fix_y[i] + "/tst_data/ktsty.txt")
    kold_array[i] *= np.loadtxt("../" + _sys_1_method_2_fix_y[i] + "/kappa_data/kappa_avg.txt")[-1]
    kold_se_array[i] = 10*kold_array[i] * np.loadtxt("../" + _sys_1_method_2_fix_y[i] + "/kappa_data/kappa_se.txt")[-1]
    knew_array[i] = np.loadtxt("../" + _sys_1_method_3_fix_y[i] + "/tst_data/ktsty.txt")
    knew_array[i] *= np.loadtxt("../" + _sys_1_method_3_fix_y[i] + "/kappa_data/kappa_avg.txt")[-1]
    knew_se_array[i] = 10*knew_array[i] * np.loadtxt("../" + _sys_1_method_3_fix_y[i] + "/kappa_data/kappa_se.txt")[-1]
for i in range(5,9):
    kold_array[i] = np.loadtxt("../" + _sys_1_method_2_fix_s[i] + "/tst_data/ktsts.txt")
    kold_array[i] *= np.loadtxt("../" + _sys_1_method_2_fix_s[i] + "/kappa_data/kappa_avg.txt")[-1]
    kold_se_array[i] = 10*kold_array[i] * np.loadtxt("../" + _sys_1_method_2_fix_s[i] + "/kappa_data/kappa_se.txt")[-1]
    knew_array[i] = np.loadtxt("../" + _sys_1_method_3_fix_s[i] + "/tst_data/ktsts.txt")
    knew_array[i] *= np.loadtxt("../" + _sys_1_method_3_fix_s[i] + "/kappa_data/kappa_avg.txt")[-1]
    knew_se_array[i] = 10*knew_array[i] * np.loadtxt("../" + _sys_1_method_3_fix_s[i] + "/kappa_data/kappa_se.txt")[-1]

#for i in range(9):
#    kold_array[i] = np.loadtxt("../" + _sys_1_method_3_fix_y[i] + "/tst_data/ktsty.txt")
#    #kold_array[i] *= np.loadtxt("../" + _sys_1_method_3_fix_y[i] + "/kappa_data/kappa_avg.txt")[-1]
#    kold_se_array[i] = 0*kold_array[i] * np.loadtxt("../" + _sys_1_method_3_fix_y[i] + "/kappa_data/kappa_se.txt")[-1]
#    knew_array[i] = np.loadtxt("../" + _sys_1_method_3_fix_s[i] + "/tst_data/ktsts.txt")
#    #knew_array[i] *= np.loadtxt("../" + _sys_1_method_3_fix_s[i] + "/kappa_data/kappa_avg.txt")[-1]
#    knew_se_array[i] =0* knew_array[i] * np.loadtxt("../" + _sys_1_method_3_fix_s[i] + "/kappa_data/kappa_se.txt")[-1]


gridspec_kw={'left':None,'bottom':None,'right':None,'top':None,'wspace':0.2,'hspace':0.2}
fig_kw={'figsize':(5.0,4.0),'dpi':150.0,'facecolor':"white",'edgecolor':"white",'linewidth':1}
fig, ((ax1)) = plt.subplots(1,1,sharex=False, sharey=False, gridspec_kw=gridspec_kw, **fig_kw)
ax1.tick_params(axis='both', which='major', direction='in', labelsize = 12, size = 4, width = 1.5)
ax1.spines['top'].set(linewidth=3)
ax1.spines['right'].set(linewidth=3)
ax1.spines['bottom'].set(linewidth=3)
ax1.spines['left'].set(linewidth=3)
ax1.set_xlabel(r"log(β$K_0$)", fontsize = 15)
ax1.set_ylabel(r"log($k_{\mathrm{ET}}$)", fontsize = 15)
ax1.errorbar(np.log10(beta * K0_array), np.log10(kGR_array), kGR_array*0, fmt='o-', markersize=3, linewidth=1, color='r', label=r'$k_\mathrm{GR}$')
ax1.errorbar(np.log10(beta * K0_array), np.log10(kBO_array), kBO_se_array / (kBO_array * np.log(10)), fmt='o-', markersize=3, linewidth=1, color='b', label=r'$k_\mathrm{BO}$')
ax1.errorbar(np.log10(beta * K0_array), np.log10(kIF_array), kIF_se_array, fmt='^-', markersize=3, linewidth=1, color='k', label=r'$k_\mathrm{IF}$')
ax1.errorbar(np.log10(beta * K0_array), np.log10(kold_array), kold_se_array / (kold_array * np.log(10)), fmt='s-', markersize=3, linewidth=1, color='g', label=r'$k_\mathrm{original}$')
ax1.errorbar(np.log10(beta * K0_array), np.log10(knew_array), knew_se_array / (knew_array * np.log(10)), fmt='*-', markersize=3, linewidth=1, color='gold', label=r'$k_\mathrm{new}$')
ax1.legend(loc='upper left', fontsize=9, frameon=False)
plt.show()

plt.savefig("fig1.png")

