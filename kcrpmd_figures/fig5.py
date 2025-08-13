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
_sys_3_method_1_hw_p1 = sorted([k for k in calc_dirs if k.startswith("_sys_3_method_1") and "_hw_1" in k], key=lambda s: float(s.split('_')[8]))
_sys_3_method_1_hw_n1 = sorted([k for k in calc_dirs if k.startswith("_sys_3_method_1") and "_hw_-1" in k], key=lambda s: float(s.split('_')[8]))
_sys_3_method_3_hw_p1 = sorted([k for k in calc_dirs if k.startswith("_sys_3_method_3") and "_hw_1" in k], key=lambda s: float(s.split('_')[8]))
_sys_3_method_3_hw_n1 = sorted([k for k in calc_dirs if k.startswith("_sys_3_method_3") and "_hw_-1" in k], key=lambda s: float(s.split('_')[8]))

with open("../" + _sys_3_method_1_hw_p1[0] + "/_model_params.txt") as f:
    model_params = eval(f.read())

with open("../" + _sys_3_method_1_hw_p1[0] + "/_control_params_dynamics.txt") as f:
    control_params = eval(f.read())

# ======= Pull in all the rate data =======
beta = units.hartree / (units.boltzmann * control_params["Temperature"])

Phw_p1 = np.loadtxt("../" + _sys_3_method_1_hw_p1[0] + "/tst_data/ktsts.txt")
ktst0 = np.loadtxt("../" + _sys_3_method_1_hw_p1[0] + "/tst_data/ktsts.txt")
kappa0_avg = np.loadtxt("../" + _sys_3_method_1[0] + "/kappa_data/kappa_avg.txt")[-1]
kappa0_se = np.loadtxt("../" + _sys_3_method_1[0] + "/kappa_data/kappa_se.txt")[-1]
kBO0 = ktst0 * kappa0_avg
kBO0_se = ktst0 * kappa0_se

leps_arr = np.array([key.split('_')[8] for key in _sys_3_method_1_hw_p1], dtype=float)[:]
K0_arr = np.array([key.split('_')[8] for key in _sys_3_method_1], dtype=float)[1:]
kGR_arr = np.zeros(K0_arr.shape)
kBO_arr = np.zeros(K0_arr.shape)
kBO_se_arr = np.zeros(K0_arr.shape)
kIF_arr = np.zeros(K0_arr.shape)
kIF_se_arr = np.zeros(K0_arr.shape)
knew_arr = np.zeros(K0_arr.shape)
knew_se_arr = np.zeros(K0_arr.shape)

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
plt.savefig('fig4.png')



