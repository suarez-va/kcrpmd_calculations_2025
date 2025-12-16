"""
Plotting script to visualize dynamical transmission trajectories from Part 3

"""

import os
import h5py
import matplotlib.pyplot as plt
import numpy as np

# Get current directory name (last component of the path)
#dir_two_up = os.path.basename(os.getcwd())
dir_two_up = os.path.basename(os.path.abspath(os.path.join(os.getcwd(), "..", "..")))
if dir_two_up.startswith("_sys_"):
    idx1 = dir_two_up.find("_sys_")
    idx2 = dir_two_up.find("_method_")
    idx3 = dir_two_up.find("_fix_")
    sys = int(dir_two_up[idx1 + len("_sys_")])
    method = int(dir_two_up[idx2 + len("_method_")])
    fix = dir_two_up[idx3 + len("_fix_")]
else:
    print("Must be in '_sys_.../libra_data/_itraj_...' directory.")
    exit()

if method == 1:
    with h5py.File("mem_data.hdf", 'r') as f:
        time = f["time/data"][:]
        Epot = f["Epot_ave/data"][:]
        Ekin = f["Ekin_ave/data"][:]
        s = f["q/data"][:,0,0]
    E = Epot + Ekin
elif method == 2 or method == 3:
    with h5py.File("mem_data.hdf", 'r') as f:
        time = f["time/data"][:]
        Epot = f["Epot_ave/data"][:]
        Ekin = f["Ekin_ave/data"][:]
        s = f["q/data"][:,0,0]
        y = f["y_aux_var/data"][:,0]
        ekin = f["ekin_aux_var/data"][:]
    E = Epot + Ekin + ekin

gridspec_kw={'left':None,'bottom':None,'right':None,'top':None,'wspace':0.2,'hspace':0.2}
if method == 1:
    fig_kw={'figsize':(5.0,6.0),'dpi':150.0,'facecolor':"white",'edgecolor':"white",'linewidth':1}
    fig, ((ax1,ax2)) = plt.subplots(2,1,sharex=False, sharey=False, gridspec_kw=gridspec_kw, **fig_kw)
    ax1.tick_params(axis='both', which='major', direction='in', labelsize = 12, size = 4, width = 1.5)
    ax1.spines['top'].set(linewidth=3)
    ax1.spines['right'].set(linewidth=3)
    ax1.spines['bottom'].set(linewidth=3)
    ax1.spines['left'].set(linewidth=3)
    ax1.set_ylabel("s (a.u.)", fontsize = 15)
    ax1.plot(time, s, color='k', label='s', linewidth=2)
    ax2.tick_params(axis='both', which='major', direction='in', labelsize = 12, size = 4, width = 1.5)
    ax2.spines['top'].set(linewidth=3)
    ax2.spines['right'].set(linewidth=3)
    ax2.spines['bottom'].set(linewidth=3)
    ax2.spines['left'].set(linewidth=3)
    ax2.legend(loc='upper right', fontsize=9, frameon=False)
    ax2.set_xlabel(r"time (a.u.)", fontsize = 15)
    ax2.set_ylabel("Energy (a.u.)", fontsize = 15)
    ax2.plot(time, Epot, color='blue', label='Potential', linewidth=2)
    ax2.plot(time, Ekin, color='red', label='Kinetic q', linewidth=2)
    ax2.plot(time, E, color='k', label='Total', linewidth=2)

if method == 2 or method == 3:
    fig_kw={'figsize':(5.0,9.0),'dpi':150.0,'facecolor':"white",'edgecolor':"white",'linewidth':1}
    fig, ((ax1,ax2,ax3)) = plt.subplots(3,1,sharex=False, sharey=False, gridspec_kw=gridspec_kw, **fig_kw)
    ax1.tick_params(axis='both', which='major', direction='in', labelsize = 12, size = 4, width = 1.5)
    ax1.spines['top'].set(linewidth=3)
    ax1.spines['right'].set(linewidth=3)
    ax1.spines['bottom'].set(linewidth=3)
    ax1.spines['left'].set(linewidth=3)
    ax1.set_ylabel("s (a.u.)", fontsize = 15)
    ax1.plot(time, s, color='k', label='s', linewidth=2)
    ax2.tick_params(axis='both', which='major', direction='in', labelsize = 12, size = 4, width = 1.5)
    ax2.spines['top'].set(linewidth=3)
    ax2.spines['right'].set(linewidth=3)
    ax2.spines['bottom'].set(linewidth=3)
    ax2.spines['left'].set(linewidth=3)
    ax2.set_ylabel("y (a.u.)", fontsize = 15)
    ax2.axhline(y=-1.5, linestyle = '--', color='k', linewidth=1)
    ax2.axhline(y=-0.5, linestyle = '--', color='k', linewidth=1)
    ax2.axhline(y=0.5, linestyle = '--', color='k', linewidth=1)
    ax2.axhline(y=1.5, linestyle = '--', color='k', linewidth=1)
    ax2.plot(time, y, color='k', label='y', linewidth=2)
    ax3.tick_params(axis='both', which='major', direction='in', labelsize = 12, size = 4, width = 1.5)
    ax3.spines['top'].set(linewidth=3)
    ax3.spines['right'].set(linewidth=3)
    ax3.spines['bottom'].set(linewidth=3)
    ax3.spines['left'].set(linewidth=3)
    ax3.legend(loc='upper right', fontsize=9, frameon=False)
    ax3.set_xlabel(r"time (a.u.)", fontsize = 15)
    ax3.set_ylabel("Energy (a.u.)", fontsize = 15)
    ax3.plot(time, Epot, color='blue', label='Potential', linewidth=2)
    ax3.plot(time, Ekin, color='red', label='Kinetic q', linewidth=2)
    ax3.plot(time, ekin, color='green', label='Kinetic y', linewidth=2)
    ax3.plot(time, E, color='k', label='Total', linewidth=2)
plt.legend()

plt.savefig('transmission_trajectory', bbox_inches='tight')
plt.show()

