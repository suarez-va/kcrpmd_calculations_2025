import sys
import os
import argparse
from liblibra_core import *
import libra_py.dynamics.tsh.compute as tsh_dynamics

# Add the parent directory to sys.path
parent_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, parent_dir)

from kcrpmd_utils.kcrpmdmodel import kcrpmd_system_bath

parser = argparse.ArgumentParser()
parser.add_argument('--nsteps', default=25000000, type=int)
parser.add_argument('--dt', default=41.34, type=float)
args = parser.parse_args()

# check current directory name
current_dir = os.path.basename(os.getcwd())
if current_dir.startswith("_sys_"):
    idx = current_dir.find("_fix_")
    fix = current_dir[idx + len("_fix_")]
else:
    print("not in correct directory, directory should start with '_sys_'.")
    exit()

# ======= READ IN KC-RPMD RECIPE, MODEL, AND INITIAL CONDITIONS =======

with open("_control_params_thermalization.txt") as f:
    control_params = eval(f.read())

with open("_model_params.txt") as f:
    model_params = eval(f.read())

with open("_init_nucl_thermalization.txt") as f:
    init_nucl = eval(f.read())

with open("_init_elec_thermalization.txt") as f:
    init_elec = eval(f.read())

control_params.update({"nsteps": args.nsteps})
control_params.update({"dt": args.dt})

rnd = Random()

res = tsh_dynamics.generic_recipe(control_params, kcrpmd_system_bath, model_params, init_elec, init_nucl, rnd)

