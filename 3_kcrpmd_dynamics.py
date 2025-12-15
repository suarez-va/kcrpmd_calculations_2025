"""
Part 3 of full calculation

This script is takes as argument:
    --itraj: index of trajectory
    --istart: starting index of configuration to be pulled from thermalization data
    --iskip: index spacing between configuration to be pulled from thermalization data separating trajectories

This script reads in the _control_params_dynamics.txt generated from 1_kcrpmd_tst.py,
then it reads in a configuration from thermalization calculation (part 2).
From the sampled configurations and the control parameters, dynamical trajectories are generated

"""

import sys
import os
import h5py
import numpy as np
import argparse

from liblibra_core import *
import util.libutil as comn
from libra_py import units
import libra_py.dynamics.tsh.compute as tsh_dynamics

# Add the parent directory to sys.path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
sys.path.insert(0, parent_dir)

from kcrpmd_utils.kcrpmdmodel import kcrpmd_system_bath

parser = argparse.ArgumentParser()
parser.add_argument('--itraj', default=1, type=int, help='transmission trajectory index')
parser.add_argument('--istart', default=999, type=int, help='thermalization starting index')
parser.add_argument('--iskip', default=9, type=int, help='thermalization skipping index')
args = parser.parse_args()

pref = F"_itraj_{args.itraj}"

# Check current directory name
current_dir = os.path.basename(os.getcwd())
if current_dir.startswith("libra_data"):
    pass
else:
    print("not in correct directory, directory should be libra_data")
    exit()

#########################################################################
# ======= READ IN KC-RPMD RECIPE, MODEL, AND INITIAL CONDITIONS ======= #
#########################################################################

with open("../_control_params_dynamics.txt") as f:
    control_params = eval(f.read())

with open("../_model_params.txt") as f:
    model_params = eval(f.read())

model_params.update({"hw": 0})
control_params.update({ "prefix":pref, "prefix2":pref })

########################################################################################
# ======= INITIALIZE TRANSMISSION TRAJECTORIES FROM THERMALIZATION CALCULATION ======= #
########################################################################################

with h5py.File("mem_data.hdf", 'r') as f:
    q = list(f["q/data"][args.itraj * args.iskip + args.istart, 0, :])
    if "use_kcrpmd" in control_params:
        y = f["y_aux_var/data"][args.itraj * args.iskip + args.istart, 0]

beta = units.hartree / (units.boltzmann * control_params["Temperature"])
mass = [model_params["ms"]] + model_params["Mj"] + [model_params["mq"]]
p = [np.random.normal(scale = np.sqrt(mass[i] / beta)) for i in range(len(mass))]

init_nucl = {"q":q, "p":p, "mass":mass, "force_constant":[0] * len(q), "init_type":0,
             "ntraj":control_params["ntraj"], "ndof": len(q)}

if "use_kcrpmd" in control_params:
    my = control_params["kcrpmd_my"]
    if my < 10.:
        my = 10.
    control_params["kcrpmd_gamma"] = np.sqrt(control_params["kcrpmd_my"] / my) * control_params["kcrpmd_gamma"]
    py = np.random.normal(scale = np.sqrt(my / beta))

    init_elec = {"init_type":0, "nstates":control_params["nstates"], "istate":0,
                  "rep":1, "ntraj":control_params["ntraj"],
                  "ndia":control_params["nstates"], "nadi":control_params["nstates"],
                  "y_aux_var":[y], "p_aux_var":[py], "m_aux_var":[my]}
else:
    init_elec = {"init_type":0, "nstates":control_params["nstates"], "istate":0,
                  "rep":1, "ntraj":control_params["ntraj"],
                  "ndia":control_params["nstates"], "nadi":control_params["nstates"]}

rnd = Random()

res = tsh_dynamics.generic_recipe(control_params, kcrpmd_system_bath, model_params, init_elec, init_nucl, rnd)

