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

from plot_utils import set_style

# Add the parent directory to sys.path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir)

from kcrpmd_utils.kcrpmdtst import KcrpmdTst

set_style()

inf = 10.
hnad = 1.
x_ar = np.linspace(-1.51, 1.51, 1000)
y_ar = (inf * (1 - np.heaviside(x_ar + 1.5, 0.5) + np.heaviside(x_ar - 1.5, 0.5))
        + hnad * (np.heaviside(x_ar + 0.5, 0.5) - np.heaviside(x_ar - 0.5, 0.5))) 

fig, ax = plt.subplots()
ax.plot(x_ar, y_ar, color='k', label='data1')
ax.set_xlim(-1.6,1.6)
ax.set_ylim(-0.1,1.5)
ax.set_xticks([-1.5, -0.5, 0.5, 1.5])
ax.set_yticks([])
ax.set_xlabel("y")
ax.set_ylabel(r"$F^{\mathrm{KC}}(y)$")
ax.set_title("")
plt.show()

