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

from kcrpmd_utils.kcrpmdtst import KcrpmdTst
from kcrpmd_utils.kcrpmdmodel import gen_kcrpmd_bath_params, get_ABC, kcrpmd_system_bath

parser = argparse.ArgumentParser()
parser.add_argument('--sys', default=1, type=int, help='KCRPMD system type A (1), B (2) or C (3)')
parser.add_argument('--fix', default='s', type=str, help='fix y or s')
parser.add_argument('--method', default=3, type=int, help='Adiabatic (1), Original KC-RPMD (2), New KC-RPMD (3)')
parser.add_argument('--a', default=0.1, type=float)
parser.add_argument('--nsteps', default=25000000, type=int)
parser.add_argument('--dt', default=41.34, type=float)
parser.add_argument('--K0', default=2.85e-3, type=float)
parser.add_argument('--leps', default=-1.43e-2, type=float)
parser.add_argument('--hw', default=0, type=int, help='left side (-1), right side (1), no hard wall (0)')
args = parser.parse_args()

if (args.sys != 1 and args.sys != 2 and args.sys != 3): print("Invalid System Type!"); exit() 
if (args.fix != 's' and args.fix != 'y'): print("Invalid Reaction Coordinate!"); exit() 
if (args.fix == 'y' and args.method == 1): print("No y Coordinate For Adiabatic Method!"); exit() 
if (args.method != 1 and args.method != 2 and args.method != 3): print("Invalid Method!"); exit() 

if (args.sys == 1 or args.sys == 2):
    if (args.method == 2):
        pref = F"_sys_{args.sys}_method_{args.method}_a_{args.a}_fix_{args.fix}_K0_{args.K0:.2e}"
    else:
        pref = F"_sys_{args.sys}_method_{args.method}_fix_{args.fix}_K0_{args.K0:.2e}"
else:
    if (args.method == 2):
        pref = F"_sys_{args.sys}_method_{args.method}_a_{args.a}_fix_{args.fix}_leps_{args.leps:.2e}_hw_{args.hw}"
    else:
        pref = F"_sys_{args.sys}_method_{args.method}_fix_{args.fix}_leps_{args.leps:.2e}_hw_{args.hw}"

os.makedirs(pref, exist_ok=True)

# ======= TST code to evaluate eta, gamma, and mass of auxiliary variable =======
T = 300 # Temperature in K
beta = units.hartree / (units.boltzmann * T)
a = args.a
b = 1000.0
c = 0.5
d = 3.0

ms = 1836.0
ws = 2.28e-3
s0 = -2.4 
s1 = 2.4 
eps = 0.0

wj, cj, mj = gen_kcrpmd_bath_params({"M":1836.0, "wc":2.28e-3, "gam":4.18608, "f":12})

K0 = args.K0 
bq = 0.0 if args.sys==1 else 3.0
mq = 5.0e4 
wq = 5.0e-4
Dq = 1.0e-4
(q0, Ea, leps) = (2.1, 6.65e-3, args.leps)
(Aq, Bq, Cq) = get_ABC(q0, leps, Ea)
qhw = 1.0
khw = 1.0e5

if args.sys == 1:
    Kq = lambda q: np.full_like(q, K0)
    Vq = lambda q: 0.5 * mq * wq**2 * q**2
elif args.sys == 2:
    Kq = lambda q: K0 * np.exp(-bq * q)
    Vq = lambda q: np.piecewise(q, [q >= 0., q < 0.],
                                [lambda q: 0.5 * mq * wq**2 * q**2,
                                 lambda q: Dq * (1 - np.exp(-np.sqrt(0.5 * mq * wq**2 / Dq) * q))**2])
elif args.sys == 3:
    Kq = lambda q: K0 * np.exp(-bq * q)
    Vq = lambda q: np.piecewise(q, [q >= 0., q < 0.],
                                [lambda q: Aq * q**4 + Bq * q**3 + Cq * q**2,
                                 lambda q: Dq * (1 - np.exp(-np.sqrt(Cq / Dq) * q))**2])

kcrpmd_tst = KcrpmdTst(beta, a, b, c, d, 1., ms, ws, s0, s1, eps, Kq, Vq)
ydag = kcrpmd_tst.ydag
sdag = kcrpmd_tst.sdag
kcrpmd_tst.set_eta_my_gammay()

if args.method == 2:
    kcrpmd_tst.eta = 2 * kcrpmd_tst.eta - np.sqrt(np.pi / kcrpmd_tst.a)
    kcrpmd_tst.a = 2 * kcrpmd_tst.a
    kcrpmd_tst.c = 0.0
    kcrpmd_tst.d = 0.0

Fg = kcrpmd_tst.Fg(); FKC = kcrpmd_tst.FKC()

if args.hw == -1:
    kcrpmd_tst.q_low = qhw - (100 / (khw * beta))**(1/6)
    Vhw = lambda q: np.piecewise(q, [q >= qhw, q < qhw], [lambda q: np.full_like(q, 0.), lambda q: khw * (q - qhw)**6])
    kcrpmd_tst.Vq = lambda q: Vq(q) + Vhw(q) 
elif args.hw == 1:
    kcrpmd_tst.q_high = qhw + (100 / (khw * beta))**(1/6)
    Vhw = lambda q: np.piecewise(q, [q >= qhw, q < qhw], [lambda q: khw * (q - qhw)**6, lambda q: np.full_like(q, 0.)])
    kcrpmd_tst.Vq = lambda q: Vq(q) + Vhw(q) 

os.makedirs(pref + "/tst_data", exist_ok=True)

#s_arr = kcrpmd_tst.s_array()
#q_arr = kcrpmd_tst.q_array()
#y_arr = kcrpmd_tst.y_array()
s_arr = np.linspace(kcrpmd_tst.s_low, kcrpmd_tst.s_high, 1000.)
q_arr = np.linspace(kcrpmd_tst.q_low, kcrpmd_tst.q_high, 1000.)
y_arr = np.linspace(kcrpmd_tst.y_low, kcrpmd_tst.y_high, 1000.)

if args.method == 1:
    Phw = np.exp(-beta * (kcrpmd_tst.Fg() - Fg)) 
    Fsq = kcrpmd_tst.Vg(s_arr, q_arr)
    Fsdag = kcrpmd_tst.Fgs(np.array([sdag]))
    Fsdagq = kcrpmd_tst.Vg(np.array([sdag]), q_arr)
    Psdagq = np.exp(-beta * (Fsdagq - Fsdag))[:,0]
    kGR = kcrpmd_tst.kGR()
    ktsts = kcrpmd_tst.kBO()
    np.savetxt(pref + "/tst_data/Phw.txt", [Phw])
    np.savetxt(pref + "/tst_data/Fsq.txt", Fsq)
    #np.savetxt(pref + "/tst_data/Fsdag.txt", Fsdag)
    #np.savetxt(pref + "/tst_data/Fsdagq.txt", Fsdagq)
    np.savetxt(pref + "/tst_data/Psdagq.txt", np.column_stack((q_arr, Psdagq)))
    np.savetxt(pref + "/tst_data/kGR.txt", [kGR])
    np.savetxt(pref + "/tst_data/ktsts.txt", [ktsts])
elif args.method == 2 or args.method == 3:
    Phw = np.exp(-beta * (kcrpmd_tst.F() - FKC)) 
    Fys = kcrpmd_tst.Fys(y_arr, s_arr)
    Fyq = kcrpmd_tst.Fyq(y_arr, q_arr)
    Fsq = kcrpmd_tst.Fsq(s_arr, q_arr)
    np.savetxt(pref + "/tst_data/Phw.txt", [Phw])
    np.savetxt(pref + "/tst_data/Fys.txt", Fys)
    np.savetxt(pref + "/tst_data/Fyq.txt", Fyq)
    np.savetxt(pref + "/tst_data/Fsq.txt", Fsq)
    if args.fix == "y":
        Fydag = kcrpmd_tst.Fy(np.array([ydag]))
        Fydags = kcrpmd_tst.Fys(np.array([ydag]), s_arr)
        Fydagq = kcrpmd_tst.Fyq(np.array([ydag]), q_arr)
        Pydags = np.exp(-beta * (Fydags - Fydag))[:,0]
        Pydagq = np.exp(-beta * (Fydagq - Fydag))[:,0]
        ktsty = kcrpmd_tst.tst_y()
        np.savetxt(pref + "/tst_data/Pydags.txt", np.column_stack((s_arr, Pydags)))
        np.savetxt(pref + "/tst_data/Pydagq.txt", np.column_stack((q_arr, Pydagq)))
        np.savetxt(pref + "/tst_data/ktsty.txt", [ktsty])
    elif args.fix == "s":
        Fsdag = kcrpmd_tst.Fs(np.array([sdag]))
        Fysdag = kcrpmd_tst.Fys(y_arr, np.array([sdag]))
        Fsdagq = kcrpmd_tst.Fsq(np.array([sdag]), q_arr)
        Pysdag = np.exp(-beta * (Fysdag - Fsdag))[0,:] 
        Psdagq = np.exp(-beta * (Fsdagq - Fsdag))[:,0] 
        ktsts = kcrpmd_tst.tst_s()
        np.savetxt(pref + "/tst_data/Pysdag.txt", np.column_stack((y_arr, Pysdag)))
        np.savetxt(pref + "/tst_data/Psdagq.txt", np.column_stack((q_arr, Psdagq)))
        np.savetxt(pref + "/tst_data/ktsts.txt", [ktsts])


# ======= SAVE LIBRA RELATED PARAMETERS  =======
nstates = 2 
ndof = 1 + len(mj) + int(args.sys != 0)
ntraj = 1

# ======= MODEL SYSTEM PARAMETERS =======
# save model parameters
_model_params = {"ms":ms, "ws":ws, "s0":s0, "s1":s1, "eps":eps,
                 "wj":wj, "cj":cj, "Mj":mj, "K0":K0,
                 "sys_type":args.sys, "mq":mq, "wq":wq, "bq":bq,
                 "Aq":Aq, "Bq":Bq, "Cq":Cq, "Dq":Dq,
                 "hard_wall":args.hw, "qhw":qhw, "khw":khw,
                 "model":1, "model0":1, "nstates": nstates}

with open(pref +  "/_model_params.txt", "w") as f:
    f.write(str(_model_params))

# ======= CHOOSE NON-ADIABATIC METHOD =======

dyn_params = {"dt":args.dt, "num_electronic_substeps":1, "nsteps":args.nsteps, "prefix":pref, "prefix2":pref,
              "hdf5_output_level":-1, "mem_output_level":3, "txt_output_level":-1,
              "use_compression":0, "compression_level":[0,0,0], "progress_frequency":0.05,
              "ntraj":ntraj, "nstates":nstates}

# General Adiabatic
def load_adiabatic(dyn_general, args):
    dyn_general.update({"rep_tdse":1}) # adiabatic representation, wfc
    dyn_general.update({"ham_update_method":1})  # recompute only diabatic Hamiltonian
    dyn_general.update({"ham_transform_method":1 }) # diabatic->adiabatic according to internal diagonalization
    dyn_general.update({"rep_force":1} ) # adiabatic
    dyn_general.update({"force_method":1} ) # state-specific  as in the TSH or adiabatic
    dyn_general.update({"time_overlap_method":1 }) # explicitly compute it from the wavefunction info
    dyn_general.update({"isNBRA":0}) # no NBRA - Hamiltonians for all trajectories are computed explicitly [ default ]
    dyn_general.update({"Temperature":args.temp}) #Temperature of the system [ default ]
    dyn_general.update({"electronic_integrator":-1}) # No propagation

# General KC-RPMD
def load_kcrpmd(dyn_general, args, kcrpmd_tst):
    dyn_general.update({"rep_tdse":0}) #diabatic representation, wfc
    dyn_general.update({"ham_update_method":1})  # recompute only diabatic Hamiltonian
    dyn_general.update({"ham_transform_method":0}) # don't do any transforms
    dyn_general.update({"rep_force":0}) # diabatic
    dyn_general.update({"force_method":4}) # KC-RPMD force 
    dyn_general.update({"time_overlap_method":1}) # explicitly compute it from the wavefunction info
    dyn_general.update({"use_kcrpmd":1}) # use it 
    dyn_general.update({"kcrpmd_a":kcrpmd_tst.a})
    dyn_general.update({"kcrpmd_b":kcrpmd_tst.b})
    dyn_general.update({"kcrpmd_c":kcrpmd_tst.c})
    dyn_general.update({"kcrpmd_d":kcrpmd_tst.d})
    dyn_general.update({"kcrpmd_eta":kcrpmd_tst.eta})
    dyn_general.update({"kcrpmd_gamma":kcrpmd_tst.gammay})
    dyn_general.update({"kcrpmd_gammaKP":0.0})
    dyn_general.update({"kcrpmd_my":kcrpmd_tst.my})
    dyn_general.update({"isNBRA":0}) # no NBRA - Hamiltonians for all trajectories are computed explicitly [ default ]
    dyn_general.update({"Temperature":args.temp}) #Temperature of the system
    dyn_general.update({"electronic_integrator":-1}) # No propagation

def load_langevin(dyn_general, _ndof):
    dyn_general.update({"ensemble":1}) #NVT
    dyn_general.update({"thermostat_params":{ "thermostat_type":"Langevin"
                                             , "Temperature":dyn_general["Temperature"]
                                             , "nu_therm":0.003 }})
    dyn_general.update({"thermostat_dofs":list(range(_ndof))})

if args.method == 1:
    load_adiabatic(dyn_params, args)
    dyn_params.update({"properties_to_save":["timestep","time","q","p","f","Epot_ave","Ekin_ave","Etot_ave"]})
elif (args.method == 2 or args.method == 3):
    load_kcrpmd(dyn_params, args, kcrpmd_tst)
    dyn_params.update({"properties_to_save":["timestep","time","q","p","f","Epot_ave","Ekin_ave","Etot_ave", "y_aux_var","p_aux_var", "f_aux_var", "ekin_aux_var"]})

# save control parameters for dynamics as copy before setting thermostat and constraints
_control_params_dynamics = dyn_params.copy()
with open(pref +  "/_control_params_dynamics.txt", "w") as f:
    f.write(str(_control_params_dynamics))

# set thermostat for thermalization
load_langevin(dyn_params, ndof)
if args.fix == 's':
    dyn_params.update({"kcrpmd_gamma":0.003})
    dyn_params.update({"kcrpmd_gammaKP":0.003})
    dyn_params.update({"constrained_dofs":[0]})
    dyn_params.update({"quantum_dofs":[]})
    dyn_params["thermostat_dofs"].pop(0)

# save control parameters for thermalization as copy
_control_params_thermalization = dyn_params.copy()
with open(pref +  "/_control_params_thermalization.txt", "w") as f:
    f.write(str(_control_params_thermalization))


