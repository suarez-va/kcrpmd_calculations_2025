import sys
import cmath
import math
import os
import h5py
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata

from liblibra_core import *
import util.libutil as comn
from libra_py import units
from libra_py import data_conv
import libra_py.dynamics.tsh.compute as tsh_dynamics
import libra_py.data_savers as data_savers

from plot_utils import set_style, add_hbar, add_abar

parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir)

from kcrpmd_utils.kcrpmdtst import KcrpmdTst
from kcrpmd_utils.kcrpmdmodel import gen_kcrpmd_bath_params, get_ABC, kcrpmd_system_bath

set_style()

# Numerically computing system C bias double well parameters Aq, Bq, and Cq, from q0, leps, and Ea
(A, B, C) = get_ABC(q0=1.0, leps=-7, Ea=3)
x_ar = np.linspace(-0.5, 1.5, 1000)
y_ar = A * x_ar**4 + B * x_ar**3 + C * x_ar**2

fig, ax = plt.subplots()
ax.plot(x_ar, y_ar, color='blue')
ax.axhline(y=0, color='k', linestyle='--', linewidth=1)
ax.set_xlim(-0.3,1.3)
ax.set_ylim(-8.0,5.0)
ax.set_xticks([0.0, 1.0])
ax.set_yticks([0])
ax.set_xticklabels([r'$0$', r'$q_0$'])
ax.set_xlabel("q")
ax.set_ylabel(r"$V_{\mathrm{DA}}(q)$")
ax.set_title("")

ax = add_abar(ax, (1.0, -7), (1.0, 0.0), label=r'$ϵ$', mutation_scale=10, label_offset=0.05, linewidth=1.5, fontsize=15)
ax = add_abar(ax, (0.39, 0.0), (0.39, 3.0), label=r'$E_\mathrm{a}$', mutation_scale=10, label_offset=0.075, linewidth=1.5, fontsize=15)

plt.tight_layout()
plt.subplots_adjust(left=0.14, right=0.98, top=0.98, bottom=0.20)
plt.savefig('fig1.png')
#plt.show()

