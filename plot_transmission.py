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

# Get current directory name (last component of the path)
current_dir = os.path.basename(os.getcwd())
if current_dir.startswith("_sys_"):
    idx1 = current_dir.find("_sys_")
    idx2 = current_dir.find("_method_")
    idx3 = current_dir.find("_fix_")
    sys = int(current_dir[idx1 + len("_sys_")])
    method = int(current_dir[idx2 + len("_method_")])
    fix = current_dir[idx3 + len("_fix_")]
else:
    print("Does not start with '_sys_'.")
    exit()

time = np.loadtxt("kappa_data/time.txt")
pos_data = np.loadtxt("kappa_data/pos.txt")
mom_data = np.loadtxt("kappa_data/mom.txt")
kappa_avg = np.loadtxt("kappa_data/kappa_avg.txt")
kappa_se = np.loadtxt("kappa_data/kappa_se.txt")

icutoff = 0
itraj = 20
plt.plot(time[icutoff:], pos_data[itraj,icutoff:])
plt.plot(time[icutoff:], pos_data[0,icutoff:])
plt.plot(time[icutoff:], pos_data[333,icutoff:])
plt.plot(time[icutoff:], pos_data[666,icutoff:])
plt.plot(time[icutoff:], pos_data[999,icutoff:])
plt.show()

gridspec_kw={'left':None,'bottom':None,'right':None,'top':None,'wspace':1.2,'hspace':0.2}
fig_kw={'figsize':(5.0,3.0),'dpi':150.0,'facecolor':"white",'edgecolor':"white",'linewidth':1}
fig, ((ax1)) = plt.subplots(1,1,sharex=False, sharey=False, gridspec_kw=gridspec_kw, **fig_kw)
ax1.tick_params(axis='both', which='major', direction='in', labelsize = 12, size = 4, width = 1.5)
ax1.set_ylabel(r"$\kappa$", fontsize = 15)
ax1.spines['top'].set(linewidth=3)
ax1.spines['right'].set(linewidth=3)
ax1.spines['bottom'].set(linewidth=3)
ax1.spines['left'].set(linewidth=3)
ax1.legend(loc='upper left', fontsize=9, frameon=False)
ax1.set_xlabel(r"$\text{q}_{\mathrm{DA}}$ coordinate", fontsize = 15)
ax1.plot(time, kappa_avg, color='k', label='exact', linewidth=2)

plt.show()

exit()


#import sys
#import cmath
#import math
#import os
#import h5py
#import matplotlib.pyplot as plt   # plots
#import numpy as np
#from scipy.interpolate import griddata
#import argparse
#
#from liblibra_core import *
#import util.libutil as comn
#from libra_py import units
#from libra_py import data_conv
#import libra_py.dynamics.tsh.compute as tsh_dynamics
#import libra_py.data_savers as data_savers
#
#method = 1
##method = 3
#if method == 1:
#    with h5py.File("mem_data.hdf", 'r') as f:
#        time = f["time/data"][:]
#        Etot = f["Etot_ave/data"][:]
#        q = f["q/data"][:,0,:]
#    E = Etot
#elif method == 2 or method == 3:
#    with h5py.File("mem_data.hdf", 'r') as f:
#        time = f["time/data"][:]
#        Etot = f["Etot_ave/data"][:]
#        q = f["q/data"][:,0,:]
#        y = f["y_aux_var/data"][:,0]
#        ekin = f["ekin_aux_var/data"][:]
#    E = Etot + ekin
#
#print(E)
#icutoff = 0
#icutoff2 = 100
##plt.plot(time[icutoff:], E[icutoff:])
##plt.plot(time[icutoff:icutoff2], q[icutoff:icutoff2,0])
#plt.plot(time[icutoff:], q[icutoff:,0])
##plt.plot(time[icutoff:], y[icutoff:])
##plt.plot(time[icutoff:], q[icutoff:,-1])
#plt.show()
##exit()
#
