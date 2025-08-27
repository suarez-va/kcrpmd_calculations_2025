import sys
import os
import matplotlib.pyplot as plt   # plots
import numpy as np

from liblibra_core import *
from libra_py import units

from plot_utils import set_style, add_hbar, add_abar

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir)

set_style()

calc_dirs = [d for d in os.listdir('../') if d.startswith('_sys_')]
_sys_2_method_1 = sorted([k for k in calc_dirs if k.startswith("_sys_2_method_1")], key=lambda s: float(s.split('_')[8]))
_sys_2_method_3_fix_y = sorted([k for k in calc_dirs if k.startswith("_sys_2_method_3_fix_y")], key=lambda s: float(s.split('_')[8]))
_sys_2_method_3_fix_s = sorted([k for k in calc_dirs if k.startswith("_sys_2_method_3_fix_s")], key=lambda s: float(s.split('_')[8]))

with open("../" + _sys_2_method_1[0] + "/_model_params.txt") as f:
    model_params = eval(f.read())

with open("../" + _sys_2_method_1[0] + "/_control_params_dynamics.txt") as f:
    control_params = eval(f.read())

# ======= Pull in all the rate data =======
beta = units.hartree / (units.boltzmann * control_params["Temperature"])

ktst0 = np.loadtxt("../" + _sys_2_method_1[0] + "/tst_data/ktsts.txt")
kappa0_se = np.loadtxt("../" + _sys_2_method_1[0] + "/kappa_data/kappa_se.txt")[-1]
kappa0_avg = np.loadtxt("../" + _sys_2_method_1[0] + "/kappa_data/kappa_avg.txt")[-1]
kBO0_se = ktst0 * kappa0_se
kBO0 = ktst0 * kappa0_avg

K0_arr = np.array([key.split('_')[8] for key in _sys_2_method_1], dtype=float)[1:]
kGR_arr = np.zeros(K0_arr.shape)
kBO_arr = np.zeros(K0_arr.shape)
kBO_se_arr = np.zeros(K0_arr.shape)
kIF_arr = np.zeros(K0_arr.shape)
kIF_se_arr = np.zeros(K0_arr.shape)
knew_arr = np.zeros(K0_arr.shape)
knew_se_arr = np.zeros(K0_arr.shape)

for i, d in enumerate(_sys_2_method_1[1:]):
    kGR_arr[i] = np.loadtxt("../" + _sys_2_method_1[i+1] + "/tst_data/kGR.txt")
    kBO_arr[i] = np.loadtxt("../" + _sys_2_method_1[i+1] + "/tst_data/ktsts.txt")
    kBO_se_arr[i] = kBO_arr[i] * np.loadtxt("../" + _sys_2_method_1[i+1] + "/kappa_data/kappa_se.txt")[-1]
    kBO_arr[i] *= np.loadtxt("../" + _sys_2_method_1[i+1] + "/kappa_data/kappa_avg.txt")[-1]
    kIF_arr[i] = kGR_arr[i] * kBO_arr[i] / (kGR_arr[i] + kBO0)
    kIF_se_arr[i] = kIF_arr[i] * kBO_se_arr[i] / kBO_arr[i]

for i in range(4):
    knew_arr[i] = np.loadtxt("../" + _sys_2_method_3_fix_y[i] + "/tst_data/ktsty.txt")
    knew_se_arr[i] = knew_arr[i] * np.loadtxt("../" + _sys_2_method_3_fix_y[i] + "/kappa_data/kappa_se.txt")[-1]
    knew_arr[i] *= np.loadtxt("../" + _sys_2_method_3_fix_y[i] + "/kappa_data/kappa_avg.txt")[-1]
for i in range(4,10):
    knew_arr[i] = np.loadtxt("../" + _sys_2_method_3_fix_s[i] + "/tst_data/ktsts.txt")
    knew_se_arr[i] = knew_arr[i] * np.loadtxt("../" + _sys_2_method_3_fix_s[i] + "/kappa_data/kappa_se.txt")[-1]
    knew_arr[i] *= np.loadtxt("../" + _sys_2_method_3_fix_s[i] + "/kappa_data/kappa_avg.txt")[-1]

fig, ax = plt.subplots(dpi=200)
ax.errorbar(np.log10(beta * K0_arr), np.log10(kGR_arr), kGR_arr*0, fmt='o-', markersize=3, linewidth=1, color='r', label=r'$k_\mathrm{GR}$', zorder=3)
ax.errorbar(np.log10(beta * K0_arr), np.log10(knew_arr), knew_se_arr / (knew_arr * np.log(10)), fmt='v-', markersize=3, linewidth=1, color='g', label=r'$k_\mathrm{new}$', zorder=4)
ax.errorbar(np.log10(beta * K0_arr), np.log10(kBO_arr), kBO_se_arr / (kBO_arr * np.log(10)), fmt='o-', markersize=3, linewidth=1, color='b', label=r'$k_\mathrm{BO}$', zorder=2)
ax.errorbar(np.log10(beta * K0_arr), np.log10(kIF_arr), kIF_se_arr / (kIF_arr * np.log(10)), fmt='^-', markersize=3, linewidth=1, color='k', label=r'$k_\mathrm{IF}$', zorder=1)
ax.set_xlim(-2.1,1.1)
ax.set_ylim(-21.0,-9.0)
ax.set_xticks([-2, -1, 0, 1])
ax.set_yticks([-20, -18, -16, -14, -12, -10])
ax.set_xlabel(r"log(β$K_0$)", fontsize = 15)
ax.set_ylabel(r"log($k_{\mathrm{ET}}$)", fontsize = 15)
ax.set_title("")
ax.legend(ncol=2, fontsize=10, loc='upper left', handletextpad=0.2,  columnspacing=0.9)

plt.subplots_adjust(left=0.18, right=0.98, top=0.98, bottom=0.17)
plt.savefig('fig5.png')

