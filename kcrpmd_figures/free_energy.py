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

beta = units.hartree / (units.boltzmann * 300.0)
a = 0.1
b = 1000.0
c = 0.5
d = 3.0

ms = 1836.0
ws = 2.28e-3
s0 = -2.4
s1 = 2.4
eps = 0.0

K0_low = 4.0e-5
#K0_high = 4.0e-3
K0_high = 7.62e-3
bq = 0.0
mq = 5e4
wq = 5e-4

Kq_low = lambda q: K0_low * np.exp(-bq * q)
Kq_high = lambda q: K0_high * np.exp(-bq * q)
Vq = lambda q: 0.5 * mq * wq**2 * q**2

adiabatic_tst_low = KcrpmdTst(beta,a,b,c,d,1.,ms,ws,s0,s1,eps,Kq_low,Vq) 
newkcrpmd_tst_low = KcrpmdTst(beta,a,b,c,d,1.,ms,ws,s0,s1,eps,Kq_low,Vq); newkcrpmd_tst_low.set_eta_my_gammay() 
oldkcrpmd_tst_low = KcrpmdTst(beta,a,b,c,d,1.,ms,ws,s0,s1,eps,Kq_low,Vq); oldkcrpmd_tst_low.set_eta_my_gammay()
oldkcrpmd_tst_low.eta = 2 * oldkcrpmd_tst_low.eta - np.sqrt(np.pi / oldkcrpmd_tst_low.a)
oldkcrpmd_tst_low.a = 2 * oldkcrpmd_tst_low.a; oldkcrpmd_tst_low.c = 0.0; oldkcrpmd_tst_low.d = 0.0

adiabatic_tst_high = KcrpmdTst(beta,a,b,c,d,1.,ms,ws,s0,s1,eps,Kq_high,Vq) 
newkcrpmd_tst_high = KcrpmdTst(beta,a,b,c,d,1.,ms,ws,s0,s1,eps,Kq_high,Vq); newkcrpmd_tst_high.set_eta_my_gammay() 
oldkcrpmd_tst_high = KcrpmdTst(beta,a,b,c,d,1.,ms,ws,s0,s1,eps,Kq_high,Vq); oldkcrpmd_tst_high.set_eta_my_gammay()
oldkcrpmd_tst_high.eta = 2 * oldkcrpmd_tst_high.eta - np.sqrt(np.pi / oldkcrpmd_tst_high.a)
oldkcrpmd_tst_high.a = 2 * oldkcrpmd_tst_high.a; oldkcrpmd_tst_high.c = 0.0; oldkcrpmd_tst_high.d = 0.0

s_arr = np.linspace(-4.0, 4.0, 1000)
# EXTRA #
q_arr = np.linspace(-1.2, 2.5, 249)
y_arr = np.linspace(-1.6, 1.6, 752)
#print(newkcrpmd_tst_low.Fthetas(0, s_arr))
#print(newkcrpmd_tst_low.Fthetaq(0, q_arr))
print(newkcrpmd_tst_low.Fy(y_arr))
#print(newkcrpmd_tst_low.Fys(y_arr, s_arr))
# EXTRA #

print(newkcrpmd_tst_low.kGR())

Fgs_low = adiabatic_tst_low.Fgs(s_arr)
newFs_low = newkcrpmd_tst_low.Fs(s_arr)
oldFs_low = oldkcrpmd_tst_low.Fs(s_arr)

Fgs_high = adiabatic_tst_high.Fgs(s_arr)
newFs_high = newkcrpmd_tst_high.Fs(s_arr)
oldFs_high = oldkcrpmd_tst_high.Fs(s_arr)

gridspec_kw={'left':None,'bottom':None,'right':None,'top':None,'wspace':0.2,'hspace':0.2}
fig_kw={'figsize':(9.0,3.0),'dpi':150.0,'facecolor':"white",'edgecolor':"white",'linewidth':1}
fig, ((ax1),(ax2)) = plt.subplots(1,2,sharex=False, sharey=False, gridspec_kw=gridspec_kw, **fig_kw)
ax1.tick_params(axis='both', which='major', direction='in', labelsize = 12, size = 4, width = 1.5)
ax1.set_ylabel("F(s)", fontsize = 15)
ax1.spines['top'].set(linewidth=3)
ax1.spines['right'].set(linewidth=3)
ax1.spines['bottom'].set(linewidth=3)
ax1.spines['left'].set(linewidth=3)
ax1.legend(loc='upper left', fontsize=9, frameon=False)
ax2.tick_params(axis='both', which='major', direction='in', labelsize = 12, size = 4, width = 1.5)
ax2.set_ylabel("F(s)", fontsize = 15)
ax2.spines['top'].set(linewidth=3)
ax2.spines['right'].set(linewidth=3)
ax2.spines['bottom'].set(linewidth=3)
ax2.spines['left'].set(linewidth=3)
ax2.legend(loc='upper right', fontsize=9, frameon=False)
ax2.set_xlabel("s coordinate", fontsize = 15)
ax1.plot(s_arr, Fgs_low, color='k', label='adiabatic', linewidth=2)
ax1.plot(s_arr, newFs_low, color='r', label='new KC-RPMD', linewidth=2)
ax1.plot(s_arr, oldFs_low, color='b', label='old KC-RPMD', linewidth=2)
ax2.plot(s_arr, Fgs_high, color='k', label='adiabatic', linewidth=2)
ax2.plot(s_arr, newFs_high, color='r', label='new KC-RPMD', linewidth=2)
ax2.plot(s_arr, oldFs_high, color='b', label='old KC-RPMD', linewidth=2)
plt.show()


