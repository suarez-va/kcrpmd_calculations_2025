"""
Part 2 of full calculation

This script reads in the _control_params_thermalization.txt generated from 1_kcrpmd_tst.py, and runs the dynamics

"""

import os
import numpy as np
from liblibra_core import *
import libra_py.dynamics.tsh.compute as tsh_dynamics


from kcrpmd_utils.kcrpmdmodel import kcrpmd_system_bath

# ======= READ IN KC-RPMD RECIPE, MODEL, AND INITIAL CONDITIONS: THERMALIZATION =======

with open("_control_params_thermalization.txt") as f:
    control_params = eval(f.read())

with open("_model_params.txt") as f:
    model_params = eval(f.read())

with open("_init_nucl_thermalization.txt") as f:
    init_nucl = eval(f.read())

with open("_init_elec_thermalization.txt") as f:
    init_elec = eval(f.read())

rnd = Random()

res = tsh_dynamics.generic_recipe(control_params, kcrpmd_system_bath, model_params, init_elec, init_nucl, rnd)


# ======= PLOT THE RESULTS =======
icutoff = 1

import h5py
import matplotlib.pyplot as plt
from scipy.integrate import trapezoid, cumulative_trapezoid

plt.rcParams.update({
    'figure.figsize': (8.0, 4.0),
    'figure.dpi': 300,
    'figure.facecolor': 'white',
    'figure.edgecolor': 'white',
    'lines.linewidth': 2,
    'axes.linewidth': 3,
    'axes.labelsize': 15,
    'axes.titlesize': 15,
    'xtick.direction': 'in',
    'xtick.top': True,
    'ytick.direction': 'in',
    'ytick.right': True,
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'xtick.major.size': 4,
    'ytick.major.size': 4,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 11,
    'legend.frameon': False,
})

# Get current directory name (last component of the path)
current_dir = os.path.basename(os.getcwd())
parent_dir = os.path.basename(os.path.dirname(os.getcwd()))

fix = current_dir[current_dir.find("_fix_") + len("_fix_")]

if parent_dir == "adiabatic":
    fig, ((ax1)) = plt.subplots(1,1)
    with h5py.File("libra_data/mem_data.hdf", 'r') as f:
        q = f["q/data"][:,0,:]
    if fix == "s":
        Psdagq_data = np.loadtxt("tst_data/Psdagq.txt"); q_arr = Psdagq_data[:,0]; Psdagq_arr = Psdagq_data[:,1]
        Isdagq_arr = cumulative_trapezoid(Psdagq_arr, q_arr, initial=0)
        q_low = q_arr[np.where(Isdagq_arr <= 0.01)[0][-1]]; q_high = q_arr[np.where(Isdagq_arr >= 0.99)[0][0]]

        ax1.set_xlabel(r"q (a.u.)")
        ax1.set_xlim([q_low-0.75*(q_high-q_low), q_high+0.75*(q_high-q_low)])
        ax1.set_ylabel("P(q)dq (a.u.)")
        ax1.hist(q[icutoff:,-1], bins=49, density=True, color='skyblue', edgecolor='black', label='Langevin')
        ax1.plot(q_arr, Psdagq_arr, color='k', label='exact')
        ax1.legend(loc='upper left')

if parent_dir == "kcrpmd_ori" or parent_dir == "kcrpmd_new":
    plt.rcParams.update({'figure.figsize': (8.0, 6.0)})
    fig, ((ax1, ax2)) = plt.subplots(2,1)
    plt.subplots_adjust(hspace=0.3)
    with h5py.File("libra_data/mem_data.hdf", 'r') as f:
        q = f["q/data"][:,0,:]
        y = f["y_aux_var/data"][:,0]
    if fix == "s":
        Pysdag_data = np.loadtxt("tst_data/Pysdag.txt"); y_arr = Pysdag_data[:,0]; Pysdag_arr = Pysdag_data[:,1]
        Iysdag_arr = cumulative_trapezoid(Pysdag_arr, y_arr, initial=0)
        y_low = y_arr[np.where(Iysdag_arr <= 0.01)[0][-1]]; y_high = y_arr[np.where(Iysdag_arr >= 0.99)[0][0]]
        Psdagq_data = np.loadtxt("tst_data/Psdagq.txt"); q_arr = Psdagq_data[:,0]; Psdagq_arr = Psdagq_data[:,1]
        Isdagq_arr = cumulative_trapezoid(Psdagq_arr, q_arr, initial=0)
        q_low = q_arr[np.where(Isdagq_arr <= 0.01)[0][-1]]; q_high = q_arr[np.where(Isdagq_arr >= 0.99)[0][0]]

        ax1.set_xlabel(r"y (a.u.)")
        ax1.set_xlim([y_low-0.75*(y_high-y_low), y_high+0.75*(y_high-y_low)])
        ax1.set_ylabel("P(y)dy (a.u.)")
        ax1.hist(y[icutoff:], bins=49, density=True, color='skyblue', edgecolor='black', label='Langevin')
        ax1.plot(y_arr, Pysdag_arr, color='k', label='exact')
        ax1.legend(loc='upper left')

        ax2.set_xlabel(r"q (a.u.)")
        ax2.set_xlim([q_low-0.75*(q_high-q_low), q_high+0.75*(q_high-q_low)])
        ax2.set_ylabel("P(q)dq (a.u.)")   
        ax2.hist(q[icutoff:,-1], bins=49, density=True, color='skyblue', edgecolor='black', label='Langevin')
        ax2.plot(q_arr, Psdagq_arr, color='k', label='exact')
        ax2.legend(loc='upper left')

    if fix == "y":
        Pydags_data = np.loadtxt("tst_data/Pydags.txt"); s_arr = Pydags_data[:,0]; Pydags_arr = Pydags_data[:,1]
        Iydags_arr = cumulative_trapezoid(Pydags_arr, s_arr, initial=0)
        s_low = s_arr[np.where(Iydags_arr <= 0.01)[0][-1]]; s_high = s_arr[np.where(Iydags_arr >= 0.99)[0][0]]
        Pydagq_data = np.loadtxt("tst_data/Pydagq.txt"); q_arr = Pydagq_data[:,0]; Pydagq_arr = Pydagq_data[:,1]
        Iydagq_arr = cumulative_trapezoid(Pydagq_arr, q_arr, initial=0)
        q_low = q_arr[np.where(Iydagq_arr <= 0.01)[0][-1]]; q_high = q_arr[np.where(Iydagq_arr >= 0.99)[0][0]]

        ax1.set_xlabel(r"s (a.u.)")
        ax1.set_xlim([s_low-0.75*(s_high-s_low), s_high+0.75*(s_high-s_low)])
        ax1.set_ylabel("P(s)ds (a.u.)")
        ax1.hist(q[icutoff:,0], bins=49, density=True, color='skyblue', edgecolor='black', label='Langevin')
        ax1.plot(s_arr, Pydags_arr, color='k', label='exact')
        ax1.legend(loc='upper left')

        ax2.set_xlabel(r"q (a.u.)")
        ax2.set_xlim([q_low-0.75*(q_high-q_low), q_high+0.75*(q_high-q_low)])
        ax2.set_ylabel("P(q)dq (a.u.)")   
        ax2.hist(q[icutoff:,-1], bins=49, density=True, color='skyblue', edgecolor='black', label='Langevin')
        ax2.plot(q_arr, Pydagq_arr, color='k', label='exact')
        ax2.legend(loc='upper left')

plt.savefig('thermalization', bbox_inches='tight')

