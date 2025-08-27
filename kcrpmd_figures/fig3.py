import sys
import os
import matplotlib.pyplot as plt   # plots
import numpy as np

from liblibra_core import *
from libra_py import units

from plot_utils import set_style

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir)

from kcrpmd_utils.kcrpmdtst import KcrpmdTst

set_style()

T = 300.0
beta = units.hartree / (units.boltzmann * T)
a = 0.1
b = 1000.0
c = 0.5
d = 3.0

ms = 1836.0
ws = 2.28e-3
s0 = -2.4
s1 = 2.4
eps = 0.0

#K0_low = 4.0e-5
#K0_high = 4.0e-3
#K0_high = 7.62e-3
K0_low = 9.55e-05
#K0_high = 4.47e-03
K0_high = 6.46e-03
#K0_high = 9.33e-03
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
oldkcrpmda1_tst_high = KcrpmdTst(beta,1.0,b,c,d,1.,ms,ws,s0,s1,eps,Kq_high,Vq); oldkcrpmda1_tst_high.set_eta_my_gammay()
oldkcrpmda1_tst_high.eta = 2 * oldkcrpmda1_tst_high.eta - np.sqrt(np.pi / oldkcrpmda1_tst_high.a)
oldkcrpmda1_tst_high.a = 2 * oldkcrpmda1_tst_high.a; oldkcrpmda1_tst_high.c = 0.0; oldkcrpmda1_tst_high.d = 0.0

s_ar = np.linspace(-5.0, 5.0, 1000)

Fgs_low = adiabatic_tst_low.Fgs(s_ar)
newFs_low = newkcrpmd_tst_low.Fs(s_ar)
oldFs_low = oldkcrpmd_tst_low.Fs(s_ar)

Fgs_high = adiabatic_tst_high.Fgs(s_ar)
newFs_high = newkcrpmd_tst_high.Fs(s_ar)
oldFs_high = oldkcrpmd_tst_high.Fs(s_ar)
oldFsa1_high = oldkcrpmda1_tst_high.Fs(s_ar)

fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(4.0, 4.0), sharex=True, sharey=False, dpi=300)
ax1.text(-0.1, 0.92, "(a)", transform=ax1.transAxes, fontsize=12)
ax2.text(-0.1, 0.92, "(b)", transform=ax2.transAxes, fontsize=12)
ax1.plot(s_ar, Fgs_low, color='k', label='adiabatic', dashes=[4, 1], zorder=3)
ax1.plot(s_ar, newFs_low, color='g', label='new', zorder=2)
ax1.plot(s_ar, oldFs_low, color='darkorange', label='ori', zorder=1)
ax2.plot(s_ar, Fgs_high, color='k', label='adiabatic', dashes=[4, 1], zorder=4)
ax2.plot(s_ar, newFs_high, color='g', label='new', zorder=3)
ax2.plot(s_ar, oldFs_high, color='darkorange', label='ori', zorder=2)
ax2.plot(s_ar, oldFsa1_high, color='blueviolet', label='ori*', zorder=1)
ax1.set_xlim(-4.7,4.7)
ax1.set_ylim(-0.002,0.03)
ax2.set_ylim(-0.002,0.03)
ax1.set_yticks([0.0, 0.01, 0.02])
ax2.set_xticks([-4, -2, 0, 2, 4])
ax2.set_yticks([0.0, 0.01, 0.02])
ax2.set_xlabel(r"$s$")
ax1.set_ylabel(r"$\Delta F(s)$")
ax2.set_ylabel(r"$\Delta F(s)$")
ax1.set_title("")

plt.subplots_adjust(hspace=0.05, left=0.18, right=0.98, top=0.98, bottom=0.12)
plt.savefig('fig3.png')

