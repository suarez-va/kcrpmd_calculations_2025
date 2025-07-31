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

with open("../" + _sys_1_method_1[0] + "/_model_params.txt") as f:
    model_params = eval(f.read())

with open("../" + _sys_1_method_1[0] + "/_control_params.txt") as f:
    control_params = eval(f.read())

# ======= TST code to evaluate interpolated rates =======
beta = units.hartree / (units.boltzmann * control_params["Temperature"])

ms = model_params["ms"]
ws = model_params["ws"]
s0 = model_params["s0"]
s1 = model_params["s1"]
eps = model_params["eps"]

K0 = model_params["K0"]
mq = model_params["mq"]
wq = model_params["wq"]

Kq = lambda q: np.full_like(q, K0)
Vq = lambda q: 0.5 * mq * wq**2 * q**2

method1_tst = KcrpmdTst(beta, 0.1, 1000., 0.5, 3.0, 1., ms, ws, s0, s1, eps, Kq, Vq)

# ======= Pull in all the rate data =======
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

for i, d in enumerate(_sys_1_method_1[1:]):
    kGR_array[i] = np.loadtxt("../" + _sys_1_method_1[i+1] + "/tst_data/kGR.txt")
    kBO_array[i] = np.loadtxt("../" + _sys_1_method_1[i+1] + "/tst_data/ktsts.txt")
    kBO_array[i] *= np.loadtxt("../" + _sys_1_method_1[i+1] + "/kappa_data/kappa_avg.txt")[-1]
    kBO_se_array[i] = kBO_array[i] * np.loadtxt("../" + _sys_1_method_1[i+1] + "/kappa_data/kappa_se.txt")[-1]
    kIF_array[i] = kGR_array[i] * kBO_array[i] / (kGR_array[i] + kBO0)
    kIF_se_array[i] = kBO_se_array[i]

print(beta)
gridspec_kw={'left':None,'bottom':None,'right':None,'top':None,'wspace':0.2,'hspace':0.2}
fig_kw={'figsize':(5.0,3.0),'dpi':150.0,'facecolor':"white",'edgecolor':"white",'linewidth':1}
fig, ((ax1)) = plt.subplots(1,1,sharex=False, sharey=False, gridspec_kw=gridspec_kw, **fig_kw)
ax1.tick_params(axis='both', which='major', direction='in', labelsize = 12, size = 4, width = 1.5)
ax1.set_ylabel("lnk", fontsize = 15)
ax1.spines['top'].set(linewidth=3)
ax1.spines['right'].set(linewidth=3)
ax1.spines['bottom'].set(linewidth=3)
ax1.spines['left'].set(linewidth=3)
ax1.legend(loc='upper left', fontsize=9, frameon=False)
ax1.set_xlabel("lnbetaK0", fontsize = 15)
ax1.errorbar(np.log10(beta * K0_array), np.log10(kGR_array), kGR_array*0, fmt='o', markersize=3, color='r', label='kGR')
ax1.errorbar(np.log10(beta * K0_array), np.log10(kBO_array), kBO_se_array / (kBO_array * np.log(10)), fmt='o', markersize=3, color='b', label='kBO')
ax1.errorbar(np.log10(beta * K0_array), np.log10(kIF_array), kIF_se_array, fmt='^', markersize=3, color='k', label='kIF')
plt.show()


