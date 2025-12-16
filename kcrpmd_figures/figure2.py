import sys
import os
import matplotlib.pyplot as plt
import numpy as np

from liblibra_core import *
from libra_py import units

from plot_utils import set_style

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir)

from kcrpmd_utils.kcrpmdtst import KcrpmdTst

set_style()

# Setting parameters for free energy plot analysis, consistent with system A in this case.
T = 300.0
beta = units.hartree / (units.boltzmann * T)
a = 0.1
b = 1000.0
c = 1.0

ms = 1836.0
ws = 2.28e-3
s0 = -2.4
s1 = 2.4
eps = 0.0

# Setting coupling constant values for free energy comparison in the two limiting cases
K0_low = 9.55e-05
K0_high = 6.46e-03
bq = 0.0
mq = 5e4
wq = 5e-4

# Functional forms of diabatic coupling and donor-acceptor potential to be used in KC-RPMD TST code
Kq_low = lambda q: K0_low * np.exp(-bq * q)
Kq_high = lambda q: K0_high * np.exp(-bq * q)
Vq = lambda q: 0.5 * mq * wq**2 * q**2

# Initializing TST objects for adiabtic, new implementation KC-RPMD, and old implementation KC-RPMD
# at low diabatic coupling
adiabatic_tst_low = KcrpmdTst(beta,a,b,c,ms,ws,s0,s1,eps,Kq_low,Vq) 
newkcrpmd_tst_low = KcrpmdTst(beta,a,b,c,ms,ws,s0,s1,eps,Kq_low,Vq); newkcrpmd_tst_low.set_eta_my_gammay() 
oldkcrpmd_tst_low = KcrpmdTst(beta,a,b,c,ms,ws,s0,s1,eps,Kq_low,Vq); oldkcrpmd_tst_low.set_eta_my_gammay()
oldkcrpmd_tst_low.eta = 2 * oldkcrpmd_tst_low.eta - np.sqrt(np.pi / oldkcrpmd_tst_low.a)
oldkcrpmd_tst_low.a = 2 * oldkcrpmd_tst_low.a; oldkcrpmd_tst_low.c = 0.0

# Initializing TST objects for adiabtic, new implementation KC-RPMD, old implementation KC-RPMD,
# and old implementation with a<=1.0 criteria at high diabatic coupling
adiabatic_tst_high = KcrpmdTst(beta,a,b,c,ms,ws,s0,s1,eps,Kq_high,Vq) 
newkcrpmd_tst_high = KcrpmdTst(beta,a,b,c,ms,ws,s0,s1,eps,Kq_high,Vq); newkcrpmd_tst_high.set_eta_my_gammay() 
oldkcrpmd_tst_high = KcrpmdTst(beta,a,b,c,ms,ws,s0,s1,eps,Kq_high,Vq); oldkcrpmd_tst_high.set_eta_my_gammay()
oldkcrpmd_tst_high.eta = 2 * oldkcrpmd_tst_high.eta - np.sqrt(np.pi / oldkcrpmd_tst_high.a)
oldkcrpmd_tst_high.a = 2 * oldkcrpmd_tst_high.a; oldkcrpmd_tst_high.c = 0.0
oldkcrpmda1_tst_high = KcrpmdTst(beta,1.0,b,c,ms,ws,s0,s1,eps,Kq_high,Vq); oldkcrpmda1_tst_high.set_eta_my_gammay()
oldkcrpmda1_tst_high.eta = 2 * oldkcrpmda1_tst_high.eta - np.sqrt(np.pi / oldkcrpmda1_tst_high.a)
oldkcrpmda1_tst_high.a = 2 * oldkcrpmda1_tst_high.a; oldkcrpmda1_tst_high.c = 0.0

# Creating array to scan along s coordinate, then computing free energies for each TST object
s_ar = np.linspace(-5.0, 5.0, 1000)

Fgs_low = adiabatic_tst_low.Fgs(s_ar)
newFs_low = newkcrpmd_tst_low.Fs(s_ar)
oldFs_low = oldkcrpmd_tst_low.Fs(s_ar)

Fgs_high = adiabatic_tst_high.Fgs(s_ar)
newFs_high = newkcrpmd_tst_high.Fs(s_ar)
oldFs_high = oldkcrpmd_tst_high.Fs(s_ar)
oldFsa1_high = oldkcrpmda1_tst_high.Fs(s_ar)

# Final plots show thermally weighted free energy difference between KC-RPMD and the ground adiabat potential
fig, (ax1, ax2) = plt.subplots(nrows=2, ncols=1, figsize=(4.0, 4.0), sharex=True, sharey=False, dpi=300)
ax1.text(-0.1, 0.92, "(a)", transform=ax1.transAxes, fontsize=14)
ax2.text(-0.1, 0.92, "(b)", transform=ax2.transAxes, fontsize=14)
ax1.plot(s_ar, beta * (newFs_low - Fgs_low), color='g', label='$\Delta F_{\mathrm{new}}$', dashes=[4, 2], zorder=2)
ax1.plot(s_ar, beta * (oldFs_low - Fgs_low), color='darkorange', label='$\Delta F_{\mathrm{ori}}$', zorder=1)
ax2.plot(s_ar, beta * (newFs_high - Fgs_high), color='g', label='$\Delta F_{\mathrm{new}}$', dashes=[4, 2], zorder=3)
ax2.plot(s_ar, beta * (oldFs_high - Fgs_high), color='darkorange', label=r'$\Delta F_{\mathrm{ori}}$', zorder=2)
ax2.plot(s_ar, beta * (oldFsa1_high - Fgs_high), color='blueviolet', label=r'$\Delta F_{\mathrm{ori}}^{*}$', zorder=1)
ax1.set_xlim(-4.7,4.7)
ax1.set_ylim(-1.58, 2.95)
ax2.set_ylim(-1.58, 2.95)
ax1.set_yticks([-1, 0, 1, 2])
ax2.set_xticks([-4, -2, 0, 2, 4])
ax2.set_yticks([-1, 0, 1, 2])
ax2.set_xlabel(r"$s$")
ax1.set_ylabel(r"$\beta\Delta F(s)$", labelpad=1, fontsize = 17.5)
ax2.set_ylabel(r"$\beta\Delta F(s)$", labelpad=1, fontsize = 17.5)
ax1.legend(fontsize=11, loc='upper right', handletextpad=0.7)
ax2.legend(fontsize=10, loc='upper right', handletextpad=0.4)
plt.legend()

plt.subplots_adjust(hspace=0.05, left=0.17, right=0.98, top=0.98, bottom=0.13)
plt.savefig('fig2.png')
#plt.show()

