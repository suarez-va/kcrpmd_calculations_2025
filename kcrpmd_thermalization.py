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
parser.add_argument('--sys', default=1, type=int, help='KCRPMD system type A (1), B (2) or C (3)')
parser.add_argument('--fix', default='None', type=str, help='fix y or s')
parser.add_argument('--nsteps', default=1000, type=int)
parser.add_argument('--dt', default=1.0, type=float)
parser.add_argument('--a', default=0.1, type=float)
parser.add_argument('--K0', default=1e-3, type=float)
parser.add_argument('--bq', default=3.0, type=float)
parser.add_argument('--leps', default=-3e-3, type=float)
parser.add_argument('--hw', default=0, type=int, help='left side (-1), right side (1), no hard wall (0)')
args = parser.parse_args()

omega, coupl, mass = gen_kcrpmd_bath_params({"M":1836.0, "wc":2.28e-3, "gam":4.18608, "f":12})
(A, B, C) = get_ABC(2.1, args.leps, 6.65030428e-3)
model_params = {"ms":1836.0, "ws":2.28e-3, "s0":-2.4, "s1":2.4, "eps":0.0,
                 "wj":omega, "cj":coupl, "Mj":mass, "K0":args.K0,
                 "sys_type":args.sys, "mq":5e4, "wq":5e-4, "bq":args.bq,
                 "Aq":A, "Bq":B, "Cq":C, "Dq":1e-3,
                 "hard_wall":args.hw, "qhw":1.0, "khw":1e5, "model":1, "model0":1, "nstates": 2}

nstates = model_params["nstates"]
ndia = nstates
nadi = nstates
ndof = 1 + len(model_params["Mj"]) + int(model_params["sys_type"] != 0)
ntraj = 1
rnd = Random()

# ======= CHOOSE NON-ADIABATIC METHOD =======
dyn_params = {"dt":args.dt, "num_electronic_substeps":1,"nsteps":args.nsteps, "prefix":"",
              "hdf5_output_level":-1, "mem_output_level":3, "txt_output_level":-1,
              "use_compression":0, "compression_level":[0,0,0], "progress_frequency":0.05,
              "properties_to_save":["timestep","time","q","p","f","Epot_ave","Ekin_ave","Etot_ave",
                                    "y_aux_var","p_aux_var", "f_aux_var", "ekin_aux_var"],
              "ntraj":ntraj, "nstates":nstates}

# General KC-RPMD
def load_kcrpmd(dyn_general):
    dyn_general.update({"rep_tdse":0}) #diabatic representation, wfc
    dyn_general.update({"ham_update_method":1})  # recompute only diabatic Hamiltonian
    dyn_general.update( {"ham_transform_method":0 }) # don't do any transforms
    dyn_general.update({"rep_force":0} ) # diabatic
    dyn_general.update({"force_method":4} ) # KC-RPMD force 
    dyn_general.update( {"time_overlap_method":1 }) # explicitly compute it from the wavefunction info
    dyn_general.update({"use_kcrpmd":1}) # use it 
    dyn_general.update({"kcrpmd_a":0.1})
    dyn_general.update({"kcrpmd_b":1000.0})
    dyn_general.update({"kcrpmd_c":0.5})
    dyn_general.update({"kcrpmd_d":3.0})
    dyn_general.update({"kcrpmd_eta":2*np.pi})
    dyn_general.update({"kcrpmd_gamma":0.0})
    dyn_general.update({"isNBRA":0}) # no NBRA - Hamiltonians for all trajectories are computed explicitly [ default ]
    dyn_general.update({"Temperature":300.0}) #Temperature of the system [ default ]
    dyn_general.update({"electronic_integrator":-1}) # No propagation

def load_langevin(dyn_general, _ndof):
    dyn_general.update({"ensemble":1}) #NVT
    dyn_general.update({"thermostat_params":{ "thermostat_type":"Langevin"
                                             , "Temperature":dyn_general["Temperature"]
                                             , "nu_therm":0.003 }})
    dyn_general.update({"thermostat_dofs":list(range(_ndof))})

load_kcrpmd(dyn_params)
dyn_params.update({"kcrpmd_a":args.a})
load_langevin(dyn_params, ndof)


# ======= Bring in TST code to evaluate eta, gamma, and mass of auxiliary variable =======
beta = units.hartree / (units.boltzmann * dyn_params["Temperature"])
a = dyn_params["kcrpmd_a"]
b = dyn_params["kcrpmd_b"]
c = dyn_params["kcrpmd_c"]
d = dyn_params["kcrpmd_d"]
e = 1.

ms = model_params["ms"]
ws = model_params["ws"]
s0 = model_params["s0"]
s1 = model_params["s1"]
eps = model_params["eps"]

K0 = model_params["K0"]
bq = model_params["bq"]
mq = model_params["mq"]
wq = model_params["wq"]
Dq = model_params["Dq"]
Aq = model_params["Aq"]
Bq = model_params["Bq"]
Cq = model_params["Cq"]
qhw = model_params["qhw"]
khw = model_params["khw"]

if model_params["sys_type"] == 0:
    Kq = lambda q: np.full_like(q, K0)
    Vq = lambda q: np.full_like(q, 0.)
elif model_params["sys_type"] == 1:
    Kq = lambda q: K0 * np.exp(-bq * q)
    Vq = lambda q: 0.5 * mq * wq**2 * q**2
elif model_params["sys_type"] == 2:
    Kq = lambda q: K0 * np.exp(-bq * q)
    Vq = lambda q: np.piecewise(q, [q >= 0., q < 0.],
                                [lambda q: 0.5 * mq * wq**2 * q**2,
                                 lambda q: Dq * (1 - np.exp(-np.sqrt(0.5 * mq * wq**2 / Dq) * q))**2])
elif model_params["sys_type"] == 3:
    Kq = lambda q: K0 * np.exp(-bq * q)
    Vq = lambda q: np.piecewise(q, [q >= 0., q < 0.],
                                [lambda q: Aq * q**4 - Bq * q**3 + Cq * q**2,
                                 lambda q: Dq * (1 - np.exp(-np.sqrt(Cq / Dq) * q))**2])

kcrpmd_tst = KcrpmdTst(beta, a, b, c, d, e, ms, ws, s0, s1, eps, Kq, Vq)
kcrpmd_tst.q_low = -1.2
kcrpmd_tst.q_high = 3.
kcrpmd_tst.set_eta_my_gammay()

eta = kcrpmd_tst.eta
gammay = kcrpmd_tst.gammay
my = kcrpmd_tst.my
sdag = kcrpmd_tst.sdag

dyn_params.update({"kcrpmd_eta":eta})
dyn_params.update({"kcrpmd_gamma":gammay})
dyn_params.update({"kcrpmd_my":my})


# ======= SAVE CONTROL PARAMETERS COPY AND COMPUTE AND TST DATA FOR LATER =======
control_params_save = dyn_params.copy()
control_params_save.update({"ensemble":0}) #NVE
control_params_save.update({"thermostat_params":{ "thermostat_type":"Langevin"
                                             , "Temperature":dyn_params["Temperature"]
                                             , "nu_therm":0.0 }})
control_params_save.update({"thermostat_dofs":[]})
control_params_save.update({"constrained_dofs":[]})

s_arr = np.linspace(-4.0, 4.0, 250)
q_arr = np.linspace(-1.2, 2., 249)
y_arr = np.linspace(-1.6, 1.6, 752)

Fydag = kcrpmd_tst.Fy(np.array([0.]))
Fydags = kcrpmd_tst.Fys(np.array([0.]), s_arr)
Fydagq = kcrpmd_tst.Fyq(np.array([0.]), q_arr)
Fsdag = kcrpmd_tst.Fs(np.array([sdag]))
Fysdag = kcrpmd_tst.Fys(y_arr, np.array([sdag]))
Fsdagq = kcrpmd_tst.Fsq(np.array([sdag]), q_arr)

Pydags = np.exp(-beta * (Fydags - Fydag))[:,0]
Pydagq = np.exp(-beta * (Fydagq - Fydag))[:,0] 
Pysdag = np.exp(-beta * (Fysdag - Fsdag))[0,:] 
Psdagq = np.exp(-beta * (Fsdagq - Fsdag))[:,0] 

Fys = kcrpmd_tst.Fys(y_arr, s_arr)
Fyq = kcrpmd_tst.Fyq(y_arr, q_arr)
Fsq = kcrpmd_tst.Fsq(s_arr, q_arr)

ktsty = kcrpmd_tst.tst_y()
ktsts = kcrpmd_tst.tst_s()


# ======= CHOOSE INITIAL CONDITIONS =======

if args.fix == 'y':
    dyn_params.update({"kcrpmd_gamma":0.0})
elif args.fix == 's':
    dyn_params.update({"kcrpmd_gamma":0.005})
    dyn_params.update({"kcrpmd_gammaKP":0.005})
    dyn_params.update({"constrained_dofs":[0]})
    dyn_params.update({"quantum_dofs":[]})
    dyn_params["thermostat_dofs"].pop(0)

nucl_params = {"q":[sdag] + [0.0] * (ndof - 2) + [qhw], "p":[0.0] * ndof, "mass":[ms] + model_params["Mj"] + [75.],
               "force_constant":[4 * (ms) / (beta**2)] + [4 * (model_params["Mj"][0]) / (beta**2)] * (ndof - 2) + [4 * (75.) / (beta**2)],
               "init_type":1, "ntraj":ntraj, "ndof": ndof}

elec_params = {"init_type":0, "nstates":nstates, "istates":[1.,0.0], "rep":0, "ntraj":ntraj, "ndia":ndia, "nadi":nadi,
               "y_aux_var":[0.0], "p_aux_var":[0.0], "m_aux_var":[500.0]}

    
pref = F"_sys_{args.sys}_fix_{args.fix}_nsteps_{args.nsteps}_dt_{args.dt}_a_{args.a}_K0_{args.K0}_bq_{args.bq}_leps_{args.leps}_hw_{args.hw}"
dyn_params.update({ "prefix":pref, "prefix2":pref })
print(F"Computing {pref}")

res = tsh_dynamics.generic_recipe(dyn_params, kcrpmd_system_bath, model_params, elec_params, nucl_params, rnd)

with open(pref + "/_control_params.txt", "w") as f:
    f.write(str(control_params_save))

if args.fix == "y":
    np.savetxt(pref + "/Pydags.txt", np.column_stack((s_arr, Pydags)))
    np.savetxt(pref + "/Pydagq.txt", np.column_stack((q_arr, Pydagq)))
    np.savetxt(pref + "/ktsty.txt", [ktsty])
elif args.fix == "s":
    np.savetxt(pref + "/Pysdag.txt", np.column_stack((y_arr, Pysdag)))
    np.savetxt(pref + "/Psdagq.txt", np.column_stack((q_arr, Psdagq)))
    np.savetxt(pref + "/ktsts.txt", [ktsts])

np.savetxt(pref + "/Fys.txt", Fys)
np.savetxt(pref + "/Fyq.txt", Fyq)
np.savetxt(pref + "/Fsq.txt", Fsq)
