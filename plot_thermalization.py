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

from kcrpmdtst import KcrpmdTst
from kcrpmdmodel import gen_kcrpmd_bath_params, get_ABC, kcrpmd_system_bath
import json

# Get current directory name (last component of the path)
current_dir = os.path.basename(os.getcwd())
if current_dir.startswith("_sys_"):
    idx = current_dir.find("_fix_")
    if idx != -1:
        char_after_fix = current_dir[idx + len("_fix_")]
        print(f"Found character: {char_after_fix}")
    else:
        print("No '_fix_' in directory name.")
else:
    print("Does not start with '_sys_'.")

if char_after_fix == "y":
    Pydags_data = np.loadtxt("Pydags.txt")
    s_arr = Pydags_data[:,0]
    Pydags = Pydags_data[:,1]
    Pydagq_data = np.loadtxt("Pydagq.txt")
    q_arr = Pydagq_data[:,0]
    Pydagq = Pydagq_data[:,1]
    ktsty = np.loadtxt("ktsty.txt")
elif char_after_fix == "s":
    Pysdag_data = np.loadtxt("Pysdag.txt")
    y_arr = Pysdag_data[:,0]
    Pysdag = Pysdag_data[:,1]
    Psdagq_data = np.loadtxt("Psdagq.txt")
    q_arr = Psdagq_data[:,0]
    Psdagq = Psdagq_data[:,1]
    ktsts = np.loadtxt("ktsts.txt")
Fys = np.loadtxt("Fys.txt")
Fyq = np.loadtxt("Fyq.txt")
Fsq = np.loadtxt("Fsq.txt")

with h5py.File("mem_data.hdf", 'r') as f:
    time = f["time/data"][:]
    q = f["q/data"][:,0,:]
    y = f["y_aux_var/data"][:,0]

gridspec_kw={'left':None,'bottom':None,'right':None,'top':None,'wspace':0.2,'hspace':0.2}
fig_kw={'figsize':(9.0,3.0),'dpi':150.0,'facecolor':"white",'edgecolor':"white",'linewidth':1}
fig, ((ax1),(ax2)) = plt.subplots(1,2,sharex=False, sharey=False, gridspec_kw=gridspec_kw, **fig_kw)
ax1.tick_params(axis='both', which='major', direction='in', labelsize = 12, size = 4, width = 1.5)
ax1.set_ylabel("Probability", fontsize = 15)
ax1.spines['top'].set(linewidth=3)
ax1.spines['right'].set(linewidth=3)
ax1.spines['bottom'].set(linewidth=3)
ax1.spines['left'].set(linewidth=3)
ax1.legend(loc='upper left', fontsize=9, frameon=False)
ax2.tick_params(axis='both', which='major', direction='in', labelsize = 12, size = 4, width = 1.5)
ax2.set_ylabel("Probability", fontsize = 15)
ax2.spines['top'].set(linewidth=3)
ax2.spines['right'].set(linewidth=3)
ax2.spines['bottom'].set(linewidth=3)
ax2.spines['left'].set(linewidth=3)
ax2.legend(loc='upper right', fontsize=9, frameon=False)
ax2.set_xlabel(r"$\text{q}_{\mathrm{DA}}$ coordinate", fontsize = 15)
ax2.hist(q[:,-1], bins=49, density = True, color='skyblue', edgecolor='black', label='Langevin')
if char_after_fix == "y":
    ax1.set_xlabel("s coordinate", fontsize = 15)
    ax1.hist(q[:,0], bins=49, density = True, color='skyblue', edgecolor='black', label='Langevin')
    ax1.plot(s_arr, Pydags, color='k', label='exact', linewidth=2)
    ax2.plot(q_arr, Pydagq, color='k', label='exact', linewidth=2)
elif char_after_fix == "s":
    ax1.set_xlabel("y coordinate", fontsize = 15)
    ax1.hist(y, bins=49, density = True, color='skyblue', edgecolor='black', label='Langevin')
    ax1.plot(y_arr, Pysdag, color='k', label='exact', linewidth=2)
    ax2.plot(q_arr, Psdagq, color='k', label='exact', linewidth=2)
plt.show()
#plt.savefig('libra', bbox_inches='tight')

