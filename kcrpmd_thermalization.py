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
parser.add_argument('--b', default=1000.0, type=float)
parser.add_argument('--c', default=0.5, type=float)
parser.add_argument('--d', default=3.0, type=float)
parser.add_argument('--temp', default=300.0, type=float, help='Temperature in K')
parser.add_argument('--nsteps', default=25000000, type=int)
parser.add_argument('--dt', default=41.34, type=float)
parser.add_argument('--K0', default=2.85e-3, type=float)
parser.add_argument('--leps', default=-1.43e-2, type=float)
parser.add_argument('--hw', default=0, type=int, help='left side (-1), right side (1), no hard wall (0)')
args = parser.parse_args()

if (args.sys != 1 and args.sys != 2 and args.sys != 3): print("Invalid System Type!"); exit() 
if (args.fix != 'y' and args.fix != 's'): print("Invalid Reaction Coordinate!"); exit() 
if (args.fix == 'y' and args.method == 1): print("No y Coordinate For Adiabatic Method!"); exit() 
if (args.method != 1 and args.method != 2 and args.method != 3): print("Invalid Method!"); exit() 

if args.sys != 3:
    pref = F"_sys_{args.sys}_method_{args.method}_fix_{args.fix}_K0_{args.K0:.2e}"
else:
    pref = F"_sys_{args.sys}_method_{args.method}_fix_{args.fix}_leps_{args.leps:.2e}_hw_{args.hw}"
os.makedirs(pref, exist_ok=True)

# ======= SET MODEL SYSTEM PARAMETERS =======
omega, coupl, mass = gen_kcrpmd_bath_params({"M":1836.0, "wc":2.28e-3, "gam":4.18608, "f":12})
if args.sys == 1: bq = 0.0
else: bq = 3.0
if args.sys == 3: (A, B, C) = get_ABC(2.1, args.leps, 6.65030428e-3)
else: (A, B, C) = (1.041e-2, -4.065e-2, 3.622e-2)
model_params = {"ms":1836.0, "ws":2.28e-3, "s0":-2.4, "s1":2.4, "eps":0.0,
                 "wj":omega, "cj":coupl, "Mj":mass, "K0":args.K0,
                 "sys_type":args.sys, "mq":5e4, "wq":5e-4, "bq":bq,
                 "Aq":A, "Bq":B, "Cq":C, "Dq":1e-3,
                 "hard_wall":args.hw, "qhw":1.0, "khw":1e5, "model":1, "model0":1, "nstates": 2}

# ======= TST code to evaluate eta, gamma, and mass of auxiliary variable =======
beta = units.hartree / (units.boltzmann * args.temp)

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
                                [lambda q: Aq * q**4 + Bq * q**3 + Cq * q**2,
                                 lambda q: Dq * (1 - np.exp(-np.sqrt(Cq / Dq) * q))**2])

kcrpmd_tst = KcrpmdTst(beta, args.a, args.b, args.c, args.d, 1., ms, ws, s0, s1, eps, Kq, Vq)
sdag = kcrpmd_tst.sdag
kcrpmd_tst.set_eta_my_gammay()
if args.method == 2:
    kcrpmd_tst.eta = 2 * kcrpmd_tst.eta - np.sqrt(np.pi / kcrpmd_tst.a)
    kcrpmd_tst.a = 2 * kcrpmd_tst.a
    kcrpmd_tst.c = 0.0
    kcrpmd_tst.d = 0.0

Fg = kcrpmd_tst.Fg(); F = kcrpmd_tst.F()
if model_params["hard_wall"] == -1:
    kcrpmd_tst.q_low = qhw - (100 / (khw * beta))**(1/6)
    Vhw = lambda q: np.piecewise(q, [q >= qhw, q < qhw], [lambda q: np.full_like(q, 0.), lambda q: khw * (q - qhw)**6])
    kcrpmd_tst.Vq = lambda q: Vq(q) + Vhw(q) 
elif model_params["hard_wall"] == 1:
    kcrpmd_tst.q_high = qhw + (100 / (khw * beta))**(1/6)
    Vhw = lambda q: np.piecewise(q, [q >= qhw, q < qhw], [lambda q: khw * (q - qhw)**6, lambda q: np.full_like(q, 0.)])
    kcrpmd_tst.Vq = lambda q: Vq(q) + Vhw(q) 
else:
    Vhw = lambda q: np.full_like(q, 0.)
if args.method == 1:
    Fghw = kcrpmd_tst.Fg()
    Phw = np.exp(-beta * (Fghw - Fg)) 
else:
    Fhw = kcrpmd_tst.F()
    Phw = np.exp(-beta * (Fhw - F)) 

s_arr = kcrpmd_tst.s_array()
q_arr = kcrpmd_tst.q_array()
y_arr = kcrpmd_tst.y_array()

if args.method == 1:
    Fsq = kcrpmd_tst.Vg(s_arr, q_arr)
    kGR = kcrpmd_tst.kGR()
    if args.fix == "s":
        Fsdag = kcrpmd_tst.Fgs(np.array([sdag]))
        Fsdagq = kcrpmd_tst.Vg(np.array([sdag]), q_arr)
        Psdagq = np.exp(-beta * (Fsdagq - Fsdag))[:,0]
        ktsts = kcrpmd_tst.kBO()
elif args.method == 2 or args.method == 3:
    Fys = kcrpmd_tst.Fys(y_arr, s_arr)
    Fyq = kcrpmd_tst.Fyq(y_arr, q_arr)
    Fsq = kcrpmd_tst.Fsq(s_arr, q_arr)
    if args.fix == "y":
        Fydag = kcrpmd_tst.Fy(np.array([0.]))
        Fydags = kcrpmd_tst.Fys(np.array([0.]), s_arr)
        Fydagq = kcrpmd_tst.Fyq(np.array([0.]), q_arr)
        Pydags = np.exp(-beta * (Fydags - Fydag))[:,0]
        Pydagq = np.exp(-beta * (Fydagq - Fydag))[:,0]
        ktsty = kcrpmd_tst.tst_y()
    elif args.fix == "s":
        Fsdag = kcrpmd_tst.Fs(np.array([sdag]))
        Fysdag = kcrpmd_tst.Fys(y_arr, np.array([sdag]))
        Fsdagq = kcrpmd_tst.Fsq(np.array([sdag]), q_arr)
        Pysdag = np.exp(-beta * (Fysdag - Fsdag))[0,:] 
        Psdagq = np.exp(-beta * (Fsdagq - Fsdag))[:,0] 
        ktsts = kcrpmd_tst.tst_s()

os.makedirs(pref + "/tst_data", exist_ok=True)
if args.method == 1:
    np.savetxt(pref + "/tst_data/Phw.txt", [Phw])
    np.savetxt(pref + "/tst_data/Fsq.txt", Fsq)
    np.savetxt(pref + "/tst_data/kGR.txt", [kGR])
    if args.fix == "s":
        np.savetxt(pref + "/tst_data/Psdagq.txt", np.column_stack((q_arr, Psdagq)))
        np.savetxt(pref + "/tst_data/ktsts.txt", [ktsts])
elif args.method == 2 or args.method == 3:
    np.savetxt(pref + "/tst_data/Phw.txt", [Phw])
    np.savetxt(pref + "/tst_data/Fys.txt", Fys)
    np.savetxt(pref + "/tst_data/Fyq.txt", Fyq)
    np.savetxt(pref + "/tst_data/Fsq.txt", Fsq)
    if args.fix == "y":
        np.savetxt(pref + "/tst_data/Pydags.txt", np.column_stack((s_arr, Pydags)))
        np.savetxt(pref + "/tst_data/Pydagq.txt", np.column_stack((q_arr, Pydagq)))
        np.savetxt(pref + "/tst_data/ktsty.txt", [ktsty])
    elif args.fix == "s":
        np.savetxt(pref + "/tst_data/Pysdag.txt", np.column_stack((y_arr, Pysdag)))
        np.savetxt(pref + "/tst_data/Psdagq.txt", np.column_stack((q_arr, Psdagq)))
        np.savetxt(pref + "/tst_data/ktsts.txt", [ktsts])

exit()

# ======= CHOOSE NON-ADIABATIC METHOD =======
nstates = model_params["nstates"]
ndia = nstates
nadi = nstates
ndof = 1 + len(model_params["Mj"]) + int(model_params["sys_type"] != 0)
ntraj = 1
rnd = Random()

dyn_params = {"dt":args.dt, "num_electronic_substeps":1,"nsteps":args.nsteps, "prefix":pref, "prefix2":pref,
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

# save control parameters as copy before setting thermostat and constraints
control_params_save = dyn_params.copy()
with open(pref + "/_control_params.txt", "w") as f:
    f.write(str(control_params_save))

# set thermostat for thermalization
load_langevin(dyn_params, ndof)
if args.fix == 's':
    dyn_params.update({"kcrpmd_gamma":0.003})
    dyn_params.update({"kcrpmd_gammaKP":0.003})
    dyn_params.update({"constrained_dofs":[0]})
    dyn_params.update({"quantum_dofs":[]})
    dyn_params["thermostat_dofs"].pop(0)

# ======= CHOOSE INITIAL CONDITIONS =======

n_therm = 1000
if args.sys == 3 and args.hw == -1:
    Kref = K0 * np.exp(-bq * 2.1)
else:
    Kref = K0
ms_therm = [max(ms*ws**2*(n_therm*args.dt/(2*np.pi))**2, 2*args.a*beta*(n_therm*ms*ws**2*(s0-s1)*args.dt/(2*np.pi*beta*Kref))**2)]
mj_therm = [model_params["Mj"][i]*omega[i]**2*(n_therm*args.dt/(2*np.pi))**2 for i in range(len(omega))]
mq_therm = [mq*wq**2*(n_therm*args.dt/(2*np.pi))**2]
mass_therm = ms_therm + mj_therm + mq_therm
force_therm = [4 * mass_therm[i] / (beta**2) for i in range(len(mass_therm))]
my_therm = 500000.0

nucl_params = {"q":[sdag] + [0.0]*(ndof-2) + [qhw], "p":[0.0]*ndof, "mass": mass_therm,
               "force_constant":force_therm, "init_type":1, "ntraj":ntraj, "ndof": ndof}

elec_params = {"init_type":0, "nstates":nstates, "rep":1, "istate": 0, "ntraj":ntraj, "ndia":ndia, "nadi":nadi,
               "y_aux_var":[0.0], "p_aux_var":[0.0], "m_aux_var":[my_therm]}

res = tsh_dynamics.generic_recipe(dyn_params, kcrpmd_system_bath, model_params, elec_params, nucl_params, rnd)

