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

N_blocks = 10

# check current directory name
current_dir = os.path.basename(os.getcwd())
if current_dir.startswith("_sys_"):
    idx = current_dir.find("_fix_")
    fix = current_dir[idx + len("_fix_")]
else:
    print("not in correct directory, directory should start with '_sys_'.")
    exit()

itraj_dirs = sorted([d for d in os.listdir('.') if d.startswith('_itraj_')])

pos_data_list = []
mom_data_list = []
# ======= READING IN TRANSMISSION TRAJECTORY DATA =======
for i, d in enumerate(itraj_dirs):
    if fix == 's':
        with h5py.File(d + "/mem_data.hdf", 'r') as f:
            time = f["time/data"][:]
            s = f["q/data"][:, 0, 0]
            ps = f["p/data"][:, 0, 0]
        pos_data_list.append(s)
        mom_data_list.append(ps)
    elif fix == 'y':
        with h5py.File(d + "/mem_data.hdf", 'r') as f:
            time = f["time/data"][:]
            y = f["y_aux_var/data"][:, 0]
            py = f["p_aux_var/data"][:, 0]
        pos_data_list.append(y)
        mom_data_list.append(py)
    #print(time.shape)
    #print(s.shape)
pos_data = np.vstack(pos_data_list)
mom_data = np.vstack(mom_data_list)

xi_pts =  pos_data.shape[0]
t_pts = pos_data.shape[1]

Num_blocks = np.zeros((N_blocks, t_pts))
Den_blocks = np.zeros(N_blocks)

block_index = np.arange(xi_pts).astype(int)
np.random.shuffle(block_index)
block_index = block_index.reshape(N_blocks, int(xi_pts / N_blocks))

for i in range(N_blocks):
    for j in range(t_pts):
        Num = 0.
        for k in block_index[i]:
            st = pos_data[k,j]
            if st >= 0.:
                Num += mom_data[k,0]
        Num_blocks[i,j] = Num

    Den = 0.
    for k in block_index[i]:
        if mom_data[k,0] >= 0.:
            Den += mom_data[k,0]

    Den_blocks[i] = Den


A = np.mean(Num_blocks, axis = 0)
B = np.mean(Den_blocks)

sA = np.std(Num_blocks, axis = 0)
sB = np.std(Den_blocks)

kappa_avg = A / B

kappa_se = kappa_avg * np.sqrt((sA / A)**2 + (sB / B)**2) / np.sqrt(N_blocks - 1)

print(kappa_avg)
print(kappa_se)

os.makedirs("kappa_data", exist_ok=True)

np.savetxt("kappa_data/time.txt", time)
np.savetxt("kappa_data/pos.txt", pos_data)
np.savetxt("kappa_data/mom.txt", mom_data)
np.savetxt("kappa_data/kappa_avg.txt", kappa_avg)
np.savetxt("kappa_data/kappa_se.txt", kappa_se)

exit()


#print(itraj_dirs)
#print(len(itraj_dirs))
#print(os.listdir())
#print(sorted(os.listdir()))
# Loop through directories named 0, 1, 2, ...
#for dirname in sorted(os.listdir())
#
## ======= READ IN MODEL AND KC-RPMD RECIPE =======
#
#with open("_model_params.txt") as f:
#    model_params = eval(f.read())
#
#with open("_control_params.txt") as f:
#    control_params = eval(f.read())
#
#model_params.update({"hw": 0})
#control_params.update({"nsteps": args.nsteps})
#control_params.update({"dt": args.dt})
#
#rnd = Random()
#
## ======= INITIALIZE TRANSMISSION TRAJECTORIES FROM THERMALIZATION CALCULATION =======
#with h5py.File("mem_data.hdf", 'r') as f:
#    q = list(f["q/data"][args.itraj * args.iskip + args.istart, 0, :])
#    if "use_kcrpmd" in control_params:
#        y = f["y_aux_var/data"][args.itraj * args.iskip + args.istart, 0]
#
#beta = units.hartree / (units.boltzmann * control_params["Temperature"])
#mass = [model_params["ms"]] + model_params["Mj"] + [model_params["mq"]]
#p = [np.random.normal(scale = np.sqrt(mass[i] / beta)) for i in range(len(mass))]
#if fix == 's':
#    p[0] = abs(p[0])
#
#nucl_params = {"q":q, "p":p, "mass":mass, "force_constant":[0] * len(q), "init_type":0,
#               "ntraj":control_params["ntraj"], "ndof": len(q)}
#
#if "use_kcrpmd" in control_params:
#    #my = control_params["kcrpmd_my"]
#    my = 10.
#    control_params["kcrpmd_gamma"] = np.sqrt(control_params["kcrpmd_my"] / my) * control_params["kcrpmd_gamma"]
#    py = np.random.normal(scale = np.sqrt(my / beta))
#    if fix == 'y':
#        py = abs(py)
#    elec_params = {"init_type":0, "nstates":control_params["nstates"], "istates":[1.0,0.0],
#                   "rep":0, "ntraj":control_params["ntraj"],
#                   "ndia":control_params["nstates"], "nadi":control_params["nstates"],
#                   "y_aux_var":[y], "p_aux_var":[py], "m_aux_var":[my]}
#else:
#    elec_params = {"init_type":0, "nstates":control_params["nstates"], "istates":[1.0,0.0],
#                   "rep":0, "ntraj":control_params["ntraj"],
#                   "ndia":control_params["nstates"], "nadi":control_params["nstates"]}
#
#pref = F"_itraj_{args.itraj}"
#control_params.update({ "prefix":pref, "prefix2":pref })
#
#res = tsh_dynamics.generic_recipe(control_params, kcrpmd_system_bath, model_params, elec_params, nucl_params, rnd)

