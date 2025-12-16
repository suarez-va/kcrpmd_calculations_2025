import sys
import os
import h5py
import matplotlib.pyplot as plt
import numpy as np

from liblibra_core import *
from libra_py import units

from plot_utils import set_style, add_hbar, add_abar

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir)

set_style()

# Reading in directories for system 3 (system C)
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

beta = units.hartree / (units.boltzmann * control_params["Temperature"])

# From directory names, create array for little epsilon driving force parameter
leps_arr = np.array([key.split('_')[10] for key in _sys_3_method_1_K0hwp1], dtype=float)[:]

# Initialize data arrays for full golden rule rate, Born-Oppenheimer rate at K0=0 and K0=2.85e-03,
# interpolation formula, new KC-RPMD implementation rates, and standard errors for all rates.
kGR_arr = np.zeros(leps_arr.shape)
kBO0_arr = np.zeros(leps_arr.shape)
kBO0_se_arr = np.zeros(leps_arr.shape)
kBO_arr = np.zeros(leps_arr.shape)
kBO_se_arr = np.zeros(leps_arr.shape)
kIF_arr = np.zeros(leps_arr.shape)
kIF_se_arr = np.zeros(leps_arr.shape)
knew_arr = np.zeros(leps_arr.shape)
knew_se_arr = np.zeros(leps_arr.shape)

# Computing full rates for kGR, kBO, and kIF by evaluating rates of each hardwall separated well,
# then averaging based on well probabilities Phw_p1 (upper left well) Phw_n1 (lower right well)
# (adiabatic Free energy probabilities)
for i, d in enumerate(_sys_3_method_1_K0hwp1):
    Phw0_p1 = np.loadtxt("../" + _sys_3_method_1_K0hwp1[i] + "/tst_data/Phw.txt")
    Phw0_n1 = np.loadtxt("../" + _sys_3_method_1_K0hwn1[i] + "/tst_data/Phw.txt")
    ktst0_p1 = np.loadtxt("../" + _sys_3_method_1_K0hwp1[i] + "/tst_data/ktsts.txt")
    ktst0_n1 = np.loadtxt("../" + _sys_3_method_1_K0hwn1[i] + "/tst_data/ktsts.txt")
    kappa0_p1 = np.loadtxt("../" + _sys_3_method_1_K0hwp1[i] + "/kappa_data/kappa_avg.txt")[-1]
    kappa0_n1 = np.loadtxt("../" + _sys_3_method_1_K0hwn1[i] + "/kappa_data/kappa_avg.txt")[-1]
    kappa0_se_p1 = 0.0
    kappa0_se_n1 = 0.0
    # So now full kBO0 rate (K0=0) is kBO0 = KBO0_p1 + KBO0_n1
    kBO0_p1 = Phw0_p1 * ktst0_p1 * kappa0_p1
    kBO0_n1 = Phw0_n1 * ktst0_n1 * kappa0_n1
    kBO0_se_p1 = Phw0_p1 * ktst0_p1 * kappa0_se_p1
    kBO0_se_n1 = Phw0_n1 * ktst0_n1 * kappa0_se_n1

    Phw_p1 = np.loadtxt("../" + _sys_3_method_1_Khwp1[i] + "/tst_data/Phw.txt")
    Phw_n1 = np.loadtxt("../" + _sys_3_method_1_Khwn1[i] + "/tst_data/Phw.txt")
    kGR_p1 = Phw_p1 * np.loadtxt("../" + _sys_3_method_1_Khwp1[i] + "/tst_data/kGR.txt")
    kGR_n1 = Phw_n1 * np.loadtxt("../" + _sys_3_method_1_Khwn1[i] + "/tst_data/kGR.txt")
    ktst_p1 = np.loadtxt("../" + _sys_3_method_1_Khwp1[i] + "/tst_data/ktsts.txt")
    ktst_n1 = np.loadtxt("../" + _sys_3_method_1_Khwn1[i] + "/tst_data/ktsts.txt")
    kappa_p1 = np.loadtxt("../" + _sys_3_method_1_Khwp1[i] + "/kappa_data/kappa_avg.txt")[-1]
    kappa_n1 = np.loadtxt("../" + _sys_3_method_1_Khwn1[i] + "/kappa_data/kappa_avg.txt")[-1]
    kappa_se_p1 = np.loadtxt("../" + _sys_3_method_1_Khwp1[i] + "/kappa_data/kappa_se.txt")[-1]
    kappa_se_n1 = np.loadtxt("../" + _sys_3_method_1_Khwn1[i] + "/kappa_data/kappa_se.txt")[-1]
    kBO_p1 = Phw_p1 * ktst_p1 * kappa_p1
    kBO_n1 = Phw_n1 * ktst_n1 * kappa_n1
    kBO_se_p1 = Phw_p1 * ktst_p1 * kappa_se_p1
    kBO_se_n1 = Phw_n1 * ktst_n1 * kappa_se_n1
    # As discussed in the KC-RPMD paper, the interpolation formula KIF is applied to each well separately
    kIF_p1 = kGR_p1 * kBO_p1 / (kGR_p1 + kBO0_p1)
    kIF_n1 = kGR_n1 * kBO_n1 / (kGR_n1 + kBO0_n1)
    kIF_se_p1 = 0
    kIF_se_n1 = 0

    # Assigning final values for full reference rates based on left and right weighted well rates
    kBO0_arr[i] = kBO0_p1 + kBO0_n1
    kBO0_se_arr[i] = np.sqrt(kBO0_se_p1**2 + kBO0_se_n1**2)
    kGR_arr[i] = kGR_p1 + kGR_n1
    kBO_arr[i] = kBO_p1 + kBO_n1
    kBO_se_arr[i] = np.sqrt(kBO_se_p1**2 + kBO_se_n1**2)
    kBO_arr[i] = kBO_p1 + kBO_n1
    kIF_arr[i] = kIF_p1 + kIF_n1
    kIF_se_arr[i] = np.sqrt(kIF_se_p1**2 + kIF_se_n1**2)

# Computing full KC-RPMD rates for hardwall well separately, then the full rates are a sum of the separate well
# contributions weighted by the probability of being found in that well (KC-RPMD Free energy probabilities)
for i, d in enumerate(_sys_3_method_3_Khwp1):
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

fig, ax = plt.subplots(dpi=200)
ax.errorbar(beta * leps_arr, np.log10(kGR_arr), kGR_arr*0, fmt='o-', markersize=3, linewidth=1, color='r', label=r'$k_\mathrm{GR}$', zorder=3)
ax.errorbar(beta * leps_arr, np.log10(knew_arr), knew_se_arr / (knew_arr * np.log(10)), fmt='v-', markersize=3, linewidth=1, color='g', label=r'$k_\mathrm{new}$', zorder=4)
ax.errorbar(beta * leps_arr, np.log10(kBO_arr), kBO_se_arr / (kBO_arr * np.log(10)), fmt='o-', markersize=3, linewidth=1, color='b', label=r'$k_\mathrm{BO}$', zorder=2)
ax.errorbar(beta * leps_arr, np.log10(kIF_arr), kIF_se_arr / (kIF_arr * np.log(10)), fmt='^-', markersize=3, linewidth=1, color='k', label=r'$k_\mathrm{IF}$', zorder=1)
ax.set_xlim(-17,1)
ax.set_ylim(-21.2,-14.2)
ax.set_xticks([-15, -10, -5, 0])
ax.set_yticks([-21, -19, -17, -15])
ax.set_xlabel(r"$βϵ$", fontsize = 17.5)
ax.set_ylabel(r"log($k_{\mathrm{ET}}$)", fontsize = 17.5)
ax.set_title("")
ax.legend(ncol=2, fontsize=10, loc='upper left', handletextpad=0.2,  columnspacing=0.9)

plt.subplots_adjust(left=0.20, right=0.98, top=0.98, bottom=0.19)
plt.savefig('fig6.png')
#plt.show()

