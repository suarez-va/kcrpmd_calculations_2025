import sys
import cmath
import math
import os
import h5py
import matplotlib.pyplot as plt   # plots
import numpy as np
from scipy.interpolate import griddata
import argparse

from liblibra_core import *
import util.libutil as comn
from libra_py import units
from libra_py import data_conv
import libra_py.dynamics.tsh.compute as tsh_dynamics
import libra_py.data_savers as data_savers

method = 3
if method == 1:
    with h5py.File("mem_data.hdf", 'r') as f:
        time = f["time/data"][:]
        Etot = f["Etot_ave/data"][:]
        q = f["q/data"][:,0,:]
    E = Etot
elif method == 2 or method == 3:
    with h5py.File("mem_data.hdf", 'r') as f:
        time = f["time/data"][:]
        Etot = f["Etot_ave/data"][:]
        q = f["q/data"][:,0,:]
        y = f["y_aux_var/data"][:,0]
        ekin = f["ekin_aux_var/data"][:]
    E = Etot + ekin

print(E)
icutoff = 0
icutoff2 = 100
#plt.plot(time[icutoff:], E[icutoff:])
#plt.plot(time[icutoff:icutoff2], q[icutoff:icutoff2,0])
plt.plot(time[icutoff:], q[icutoff:,0])
plt.plot(time[icutoff:], y[icutoff:])
#plt.plot(time[icutoff:], q[icutoff:,-1])
plt.show()
#exit()

