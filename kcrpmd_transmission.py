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

# Add the parent directory to sys.path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir)

from kcrpmd_utils.kcrpmdmodel import kcrpmd_system_bath

parser = argparse.ArgumentParser()
parser.add_argument('--itraj', default=1, type=int, help='transmission trajectory index')
parser.add_argument('--istart', default=1000000, type=int, help='thermalization starting index')
parser.add_argument('--iskip', default=23999, type=int, help='thermalization skipping index')
parser.add_argument('--nsteps', default=25000, type=int)
parser.add_argument('--dt', default=0.08268, type=float)
args = parser.parse_args()

# check current directory name
current_dir = os.path.basename(os.getcwd())
if current_dir.startswith("_sys_"):
    idx = current_dir.find("_fix_")
    fix = current_dir[idx + len("_fix_")]
else:
    print("not in correct directory, directory should start with '_sys_'.")
    exit()

# ======= READ IN MODEL AND KC-RPMD RECIPE =======

with open("_model_params.txt") as f:
    model_params = eval(f.read())

with open("_control_params.txt") as f:
    control_params = eval(f.read())

model_params.update({"hw": 0})
control_params.update({"nsteps": args.nsteps})
control_params.update({"dt": args.dt})

control_params.update({"properties_to_save":["timestep","time","q","p","f","Epot_ave","Ekin_ave","Etot_ave","Cdia","Cadi"]})

rnd = Random()

# ======= INITIALIZE TRANSMISSION TRAJECTORIES FROM THERMALIZATION CALCULATION =======
with h5py.File("mem_data.hdf", 'r') as f:
    q = list(f["q/data"][args.itraj * args.iskip + args.istart, 0, :])
    if "use_kcrpmd" in control_params:
        y = f["y_aux_var/data"][args.itraj * args.iskip + args.istart, 0]

beta = units.hartree / (units.boltzmann * control_params["Temperature"])
mass = [model_params["ms"]] + model_params["Mj"] + [model_params["mq"]]
p = [np.random.normal(scale = np.sqrt(mass[i] / beta)) for i in range(len(mass))]
if fix == 's':
    p[0] = abs(p[0])

nucl_params = {"q":q, "p":p, "mass":mass, "force_constant":[0] * len(q), "init_type":0,
               "ntraj":control_params["ntraj"], "ndof": len(q)}

if "use_kcrpmd" in control_params:
    #my = control_params["kcrpmd_my"]
    my = 10.
    control_params["kcrpmd_gamma"] = np.sqrt(control_params["kcrpmd_my"] / my) * control_params["kcrpmd_gamma"]
    py = np.random.normal(scale = np.sqrt(my / beta))
    if fix == 'y':
        py = abs(py)
    elec_params = {"init_type":0, "nstates":control_params["nstates"], "istates":[1.0,0.0],
                   "rep":0, "ntraj":control_params["ntraj"],
                   "ndia":control_params["nstates"], "nadi":control_params["nstates"],
                   "y_aux_var":[y], "p_aux_var":[py], "m_aux_var":[my]}
else:
    elec_params = {"init_type":0, "nstates":control_params["nstates"], "istates":[1.0,0.0],
                   "rep":0, "ntraj":control_params["ntraj"],
                   "ndia":control_params["nstates"], "nadi":control_params["nstates"]}

pref = F"_itraj_{args.itraj}"
control_params.update({ "prefix":pref, "prefix2":pref })

res = tsh_dynamics.generic_recipe(control_params, kcrpmd_system_bath, model_params, elec_params, nucl_params, rnd)

