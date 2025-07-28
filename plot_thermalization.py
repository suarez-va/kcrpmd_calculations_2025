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

#import json

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


if method == 1:
    with h5py.File("mem_data.hdf", 'r') as f:
        time = f["time/data"][:]
        q = f["q/data"][:,0,:]
    Fsq = np.loadtxt("tst_data/Fsq.txt")
    if fix == "s":
        Psdagq_data = np.loadtxt("tst_data/Psdagq.txt")
        q_arr = Psdagq_data[:,0]
        Psdagq = Psdagq_data[:,1]
        ktsts = np.loadtxt("tst_data/ktsts.txt")
elif method == 2 or method == 3:
    with h5py.File("mem_data.hdf", 'r') as f:
        time = f["time/data"][:]
        q = f["q/data"][:,0,:]
        y = f["y_aux_var/data"][:,0]
    Fys = np.loadtxt("tst_data/Fys.txt")
    Fyq = np.loadtxt("tst_data/Fyq.txt")
    Fsq = np.loadtxt("tst_data/Fsq.txt")
    if fix == "y":
        Pydags_data = np.loadtxt("tst_data/Pydags.txt")
        s_arr = Pydags_data[:,0]
        Pydags = Pydags_data[:,1]
        Pydagq_data = np.loadtxt("tst_data/Pydagq.txt")
        q_arr = Pydagq_data[:,0]
        Pydagq = Pydagq_data[:,1]
        ktsty = np.loadtxt("tst_data/ktsty.txt")
    elif fix == "s":
        Pysdag_data = np.loadtxt("tst_data/Pysdag.txt")
        y_arr = Pysdag_data[:,0]
        Pysdag = Pysdag_data[:,1]
        Psdagq_data = np.loadtxt("tst_data/Psdagq.txt")
        q_arr = Psdagq_data[:,0]
        Psdagq = Psdagq_data[:,1]
        ktsts = np.loadtxt("tst_data/ktsts.txt")

#icutoff = 24999000
#print(time[icutoff:].shape)
icutoff = 1000000
#plt.plot(time[icutoff:], q[icutoff:,0])
plt.plot(time[icutoff:], y[icutoff:])
#plt.plot(time[icutoff:], q[icutoff:,-1])
plt.show()
#exit()
if method == 1: 
    gridspec_kw={'left':None,'bottom':None,'right':None,'top':None,'wspace':0.2,'hspace':0.2}
    fig_kw={'figsize':(5.0,3.0),'dpi':150.0,'facecolor':"white",'edgecolor':"white",'linewidth':1}
    fig, ((ax1)) = plt.subplots(1,1,sharex=False, sharey=False, gridspec_kw=gridspec_kw, **fig_kw)
    ax1.tick_params(axis='both', which='major', direction='in', labelsize = 12, size = 4, width = 1.5)
    ax1.set_ylabel("Probability", fontsize = 15)
    ax1.spines['top'].set(linewidth=3)
    ax1.spines['right'].set(linewidth=3)
    ax1.spines['bottom'].set(linewidth=3)
    ax1.spines['left'].set(linewidth=3)
    ax1.legend(loc='upper left', fontsize=9, frameon=False)
    ax1.set_xlabel(r"$\text{q}_{\mathrm{DA}}$ coordinate", fontsize = 15)
    ax1.hist(q[icutoff:,-1], bins=49, density = True, color='skyblue', edgecolor='black', label='Langevin')
    ax1.plot(q_arr, Psdagq, color='k', label='exact', linewidth=2)
else:
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
    #ax1.set_xlim([-0.01, 0.01])
    ax2.tick_params(axis='both', which='major', direction='in', labelsize = 12, size = 4, width = 1.5)
    ax2.set_ylabel("Probability", fontsize = 15)
    ax2.spines['top'].set(linewidth=3)
    ax2.spines['right'].set(linewidth=3)
    ax2.spines['bottom'].set(linewidth=3)
    ax2.spines['left'].set(linewidth=3)
    ax2.legend(loc='upper right', fontsize=9, frameon=False)
    ax2.set_xlabel(r"$\text{q}_{\mathrm{DA}}$ coordinate", fontsize = 15)
    ax2.hist(q[icutoff:,-1], bins=49, density = True, color='skyblue', edgecolor='black', label='Langevin')
    if fix == "y":
        ax1.set_xlabel("s coordinate", fontsize = 15)
        ax1.hist(q[icutoff:,0], bins=49, density = True, color='skyblue', edgecolor='black', label='Langevin')
        ax1.plot(s_arr, Pydags, color='k', label='exact', linewidth=2)
        ax2.plot(q_arr, Pydagq, color='k', label='exact', linewidth=2)
    elif fix == "s":
        ax1.set_xlabel("y coordinate", fontsize = 15)
        ax1.hist(y[icutoff:], bins=49, density = True, color='skyblue', edgecolor='black', label='Langevin')
        ax1.plot(y_arr, Pysdag, color='k', label='exact', linewidth=2)
        ax2.plot(q_arr, Psdagq, color='k', label='exact', linewidth=2)
plt.show()

#plt.savefig('libra', bbox_inches='tight')

