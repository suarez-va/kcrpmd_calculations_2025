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

parser = argparse.ArgumentParser()
parser.add_argument('--nsteps', default=1000, type=int)
parser.add_argument('--dt', default=1.0, type=float)
parser.add_argument('--isample', default=1000, type=int)
args = parser.parse_args()


# ======= READ IN MODEL AND KC-RPMD RECIPE =======

with open("_model_params.txt") as f:
    model_params = eval(f.read())

with open("_control_params.txt") as f:
    control_params = eval(f.read())

control_params.update({"nsteps": args.nsteps})
control_params.update({"dt": args.dt})

rnd = Random()

# ======= INITIALIZE TRANSMISSION TRAJECTORIES FROM THERMALIZATION CALCULATION =======

with h5py.File("mem_data.hdf", 'r') as f:
        q = list(f["q/data"][args.isample,0,:])
        p = list(f["p/data"][args.isample,0,:])
        y = f["y_aux_var/data"][args.isample,0]
        py = f["p_aux_var/data"][args.isample,0]

nucl_params = {"q":q, "p":p, "mass":[model_params["ms"]] + model_params["Mj"] + [model_params["ms"]],
               "force_constant":[0] * len(q),
               "init_type":0, "ntraj":control_params["ntraj"], "ndof": len(q)}

elec_params = {"init_type":0, "nstates":control_params["nstates"], "istates":[1.,0.0],
               "rep":0, "ntraj":control_params["ntraj"],
               "ndia":control_params["nstates"], "nadi":control_params["nstates"],
               "y_aux_var":[y], "p_aux_var":[py], "m_aux_var":[control_params["kcrpmd_my"]]}
    
pref = F"_nsteps_{args.nsteps}_dt_{args.dt}_isample_{args.isample}"
control_params.update({ "prefix":pref, "prefix2":pref })
print(F"Computing {pref}")

res = tsh_dynamics.generic_recipe(control_params, kcrpmd_system_bath, model_params, elec_params, nucl_params, rnd)


