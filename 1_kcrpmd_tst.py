"""
Part 1 of full calculation

This script is takes as argument:
    --sys: the system type (A, B or C, see KC-RPMD paper 2025)
    --fix: which reaction coordinate to fix and evaluate rate constants from (either y or s)
    --method: whether to run 1 = adiabatic, 2 = original KC-RPMD, 3 = new KC-RPMD
    --a: KC-RPMD gaussian restraint parameter (no larger than 0.1, large enough to converge free energy of kinked-pair formation)
    --K0: diabatic coupling constant (sys A) or prefactor (sys B, C)
    --leps: donor-acceptor coordinate driving force (sys C only)
    --hw: whether to include left side hard wall (-1), right side hard wall (1), or no hard wall (0)

    From kcrpmd_utils/kcrpmdtst.py code, free energies and KC-RPMD parameters for "eta",
mass of auxiliary variable "my", and Langevin frictional coefficient "gammay",
are computed numerically and saved to _sys_*/tst_data/ for later.

    Finally, Libra control parameters dictionaries are generated and saved to _sys_*/
thermalization calculation (part 2) and for transmission coefficient dynamic trajectories (part 3) 

"""

import os
import numpy as np
import argparse

from liblibra_core import *
from libra_py import units

from kcrpmd_utils.kcrpmdtst import KcrpmdTst
from kcrpmd_utils.kcrpmdmodel import gen_kcrpmd_bath_params, get_ABC

######################################################
# ======= ARGUMENT PARSER, SEE TOP OF SCRIPT ======= #
######################################################
parser = argparse.ArgumentParser()
parser.add_argument('--sys', default=1, type=int, help='KCRPMD system type A (1), B (2) or C (3)')
parser.add_argument('--fix', default='s', type=str, help='fix y or s')
parser.add_argument('--method', default=3, type=int, help='Adiabatic (1), Original KC-RPMD (2), New KC-RPMD (3)')
parser.add_argument('--a', default=0.1, type=float)
parser.add_argument('--K0', default=2.85e-3, type=float)
parser.add_argument('--leps', default=-1.43e-2, type=float)
parser.add_argument('--hw', default=0, type=int, help='left side (-1), right side (1), no hard wall (0)')
args = parser.parse_args()

if (args.sys != 1 and args.sys != 2 and args.sys != 3): print("Invalid System Type!"); exit() 
if (args.fix != 's' and args.fix != 'y'): print("Invalid Reaction Coordinate!"); exit() 
if (args.fix == 'y' and args.method == 1): print("No y Coordinate For Adiabatic Method!"); exit() 
if (args.method != 1 and args.method != 2 and args.method != 3): print("Invalid Method!"); exit() 

###########################################################################
# ======= CREATING WORKING DIRECTORIES FOR CALCULATIONS TO BE RUN ======= #
###########################################################################
if (args.sys == 1 or args.sys == 2):
    if (args.method == 2):
        pref = F"_sys_{args.sys}_method_{args.method}_a_{args.a}_fix_{args.fix}_K0_{args.K0:.2e}"
    else:
        pref = F"_sys_{args.sys}_method_{args.method}_fix_{args.fix}_K0_{args.K0:.2e}"
else:
    if (args.method == 1):
        pref = F"_sys_{args.sys}_method_{args.method}_fix_{args.fix}_K0_{args.K0:.2e}_leps_{args.leps:.2e}_hw_{args.hw}"
    elif (args.method == 2):
        pref = F"_sys_{args.sys}_method_{args.method}_a_{args.a}_fix_{args.fix}_leps_{args.leps:.2e}_hw_{args.hw}"
    else:
        pref = F"_sys_{args.sys}_method_{args.method}_fix_{args.fix}_leps_{args.leps:.2e}_hw_{args.hw}"

os.makedirs(pref, exist_ok=True)

###################################################
# ======= ASSIGNING PARAMETERS FOR SYSTEM ======= #
###################################################
T = 300.0 # Temperature in K
beta = units.hartree / (units.boltzmann * T)
a = args.a # KC-RPMD parameter a
b = 1000.0 # KC-RPMD parameter b
c = 1.0 # KC-RPMD parameter c

ms = 1836.0 # Mass of s
ws = 2.28e-3 # Frequency of s
s0 = -2.4 # Diabat 0 parabola minima
s1 = 2.4 # Diabat 1 parabola minima
eps = 0.0 # Diabat 0 to diabat 1 driving force

wj, cj, mj = gen_kcrpmd_bath_params({"M":1836.0, "wc":2.28e-3, "gam":4.18608, "f":12}) # Ohmic spectral density bath parameters

K0 = args.K0 # Diabatic coupling constant/prefactor
bq = 0.0 if (args.sys==1 or K0<=1e-10) else 3.0 # Diabatic coupling q coordinate exponential dependence
mq = 5.0e4 # q coordinate mass
wq = 5.0e-4 # q coordinate frequency
Dq = 1.0e-4 if args.sys==2 else 1.0e-3 # q coordinate morse potential parameter
(q0, Ea, leps) = (2.1, 6.65e-3, args.leps) # System C q coordinate double well parameters
(Aq, Bq, Cq) = get_ABC(q0, leps, Ea) # Numerically computing Aq, Bq, Cq from q0, leps, Ea
qhw = 1.0 # Hard wall potential location
khw = 1.0e5 # Hard wall potential strength

###########################################################################################
# ======= NOW CALLING ON TST CODE TO COMPUTE FREE ENERGIES AND KC-RPMD PARAMETERS ======= #
###########################################################################################

# Defining Kq and Vq functions of system for TST code
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

# Instantiation tst code object, computing eta, my, and gammay
kcrpmd_tst = KcrpmdTst(beta, a, b, c, ms, ws, s0, s1, eps, Kq, Vq)
ydag = kcrpmd_tst.ydag
sdag = kcrpmd_tst.sdag
# IMPORTANT! without the hardwall separation we use eta, my, gammay computed for the whole system without hardwall potential.
# with the hard wall (sys C), instead eta, my, and gammay are computed from the right side lower well with lower diabatic coupling.
if args.hw == 0:
    kcrpmd_tst.set_eta_my_gammay()
else:
    q_low_cp = kcrpmd_tst.q_low
    kcrpmd_tst.q_low = qhw
    kcrpmd_tst.set_eta_my_gammay()
    kcrpmd_tst.q_low = q_low_cp

# This is a cheap and dirty fix to recover original KC-RPMD using new KC-RPMD formalism.
if args.method == 2:
    kcrpmd_tst.eta = 2 * kcrpmd_tst.eta - np.sqrt(np.pi / kcrpmd_tst.a)
    kcrpmd_tst.a = 2 * kcrpmd_tst.a
    kcrpmd_tst.c = 0.0

# Computing ground state free energy and KC-RPMD free energy (fully integrated)
Fg = kcrpmd_tst.Fg(); FKC = kcrpmd_tst.F()

# Computing full fermi-golden rule rate saving to _sys_*/tst_data (no hardwall influence)
os.makedirs(pref + "/tst_data", exist_ok=True)
kGR_full = kcrpmd_tst.kGR()
np.savetxt(pref + "/tst_data/kGR_full.txt", [kGR_full])

# Now we add hardwall potential to evaluate free energies and rates of each half separately (sys C only)
if args.hw == -1:
    kcrpmd_tst.q_low = qhw - (100 / (khw * beta))**(1/6)
    Vhw = lambda q: np.piecewise(q, [q >= qhw, q < qhw], [lambda q: np.full_like(q, 0.), lambda q: khw * (q - qhw)**6])
    kcrpmd_tst.Vq = lambda q: Vq(q) + Vhw(q) 
elif args.hw == 1:
    kcrpmd_tst.q_high = qhw + (100 / (khw * beta))**(1/6)
    Vhw = lambda q: np.piecewise(q, [q >= qhw, q < qhw], [lambda q: khw * (q - qhw)**6, lambda q: np.full_like(q, 0.)])
    kcrpmd_tst.Vq = lambda q: Vq(q) + Vhw(q) 

# Coordinate arrays used in numerical integration and plotting for later
s_arr = kcrpmd_tst.s_array()
q_arr = kcrpmd_tst.q_array()
y_arr = kcrpmd_tst.y_array()

# Now we calculate absolutely everything that might be useful later
# Hardwall influence if it's turned on, sys C only
if args.method == 1:
    Phw = np.exp(-beta * (kcrpmd_tst.Fg() - Fg)) # Hardwall potential probability from 0 to 1
    Fsq = kcrpmd_tst.Vg(s_arr, q_arr) # Adiabatic free energy along s and q
    Fsdag = kcrpmd_tst.Fgs(np.array([sdag])) # Adiabatic free energy at sdagger
    Fsdagq = kcrpmd_tst.Vg(np.array([sdag]), q_arr) # Adiabatic free energy along q at sdagger
    Psdagq = np.exp(-beta * (Fsdagq - Fsdag))[:,0] # Adiabatic probability along q at sdagger
    kGR = kcrpmd_tst.kGR() # Fermi-golden rule rate
    ktsts = kcrpmd_tst.kBO() # Adiabatic TST rate along s
    np.savetxt(pref + "/tst_data/Phw.txt", [Phw])
    np.savetxt(pref + "/tst_data/Fsq.txt", Fsq)
    np.savetxt(pref + "/tst_data/Psdagq.txt", np.column_stack((q_arr, Psdagq)))
    np.savetxt(pref + "/tst_data/kGR.txt", [kGR])
    np.savetxt(pref + "/tst_data/ktsts.txt", [ktsts])
elif args.method == 2 or args.method == 3:
    Phw = np.exp(-beta * (kcrpmd_tst.F() - FKC)) # Hardwall potential probability from 0 to 1
    Fys = kcrpmd_tst.Fys(y_arr, s_arr) # KC-RPMD free energy along y and s
    Fyq = kcrpmd_tst.Fyq(y_arr, q_arr) # KC-RPMD free energy along y and q
    Fsq = kcrpmd_tst.Fsq(s_arr, q_arr) # KC-RPMD free energy along s and q
    np.savetxt(pref + "/tst_data/Phw.txt", [Phw])
    np.savetxt(pref + "/tst_data/Fys.txt", Fys)
    np.savetxt(pref + "/tst_data/Fyq.txt", Fyq)
    np.savetxt(pref + "/tst_data/Fsq.txt", Fsq)
    if args.fix == "y":
        Fydag = kcrpmd_tst.Fy(np.array([ydag])) # KC-RPMD free energy at ydagger
        np.savetxt(pref + "/tst_data/Fydag.txt", [Fydag])
        Fydags = kcrpmd_tst.Fys(np.array([ydag]), s_arr) # KC-RPMD free energy along s at ydagger
        Fydagq = kcrpmd_tst.Fyq(np.array([ydag]), q_arr) # KC-RPMD free energy along q at ydagger
        Pydags = np.exp(-beta * (Fydags - Fydag))[:,0] # KC-RPMD probability along s at ydagger
        Pydagq = np.exp(-beta * (Fydagq - Fydag))[:,0] # KC-RPMD probability along q at ydagger
        ktsty = kcrpmd_tst.tst_y() # KC-RPMD TST rate along y
        np.savetxt(pref + "/tst_data/Pydags.txt", np.column_stack((s_arr, Pydags)))
        np.savetxt(pref + "/tst_data/Pydagq.txt", np.column_stack((q_arr, Pydagq)))
        np.savetxt(pref + "/tst_data/ktsty.txt", [ktsty])
    elif args.fix == "s":
        Fsdag = kcrpmd_tst.Fs(np.array([sdag])) # KC-RPMD free energy at sdagger
        Fysdag = kcrpmd_tst.Fys(y_arr, np.array([sdag])) # KC-RPMD free energy along y at sdagger
        Fsdagq = kcrpmd_tst.Fsq(np.array([sdag]), q_arr) # KC-RPMD free energy along q at sdagger
        Pysdag = np.exp(-beta * (Fysdag - Fsdag))[0,:] # KC-RPMD probability along y at sdagger
        Psdagq = np.exp(-beta * (Fsdagq - Fsdag))[:,0] # KC-RPMD probability along q at sdagger
        ktsts = kcrpmd_tst.tst_s() # KC-RPMD TST rate along s
        np.savetxt(pref + "/tst_data/Pysdag.txt", np.column_stack((y_arr, Pysdag)))
        np.savetxt(pref + "/tst_data/Psdagq.txt", np.column_stack((q_arr, Psdagq)))
        np.savetxt(pref + "/tst_data/ktsts.txt", [ktsts])

######################################################################################################################
# ======= NOW CREATING AND SAVING LIBRA CONTROL PARAMETERS FOR THERMALIZATION (PART 2) AND DYNAMICS (PART 3) ======= #
######################################################################################################################
nstates = 2
ndia = 2
nadi = 2
ndof = 1 + len(mj) + int(args.sys != 0)
ntraj = 1

# Save model parameters
_model_params = {"ms":ms, "ws":ws, "s0":s0, "s1":s1, "eps":eps,
                 "wj":wj, "cj":cj, "Mj":mj, "K0":K0,
                 "sys_type":args.sys, "mq":mq, "wq":wq, "bq":bq,
                 "Aq":Aq, "Bq":Bq, "Cq":Cq, "Dq":Dq,
                 "hard_wall":args.hw, "qhw":qhw, "khw":khw,
                 "model":1, "model0":1, "nstates": nstates}

with open(pref +  "/_model_params.txt", "w") as f:
    f.write(str(_model_params))

# Default parameters for thermalization calculation (dt, nsteps, and nprint will change for dynamics)
dyn_params = {"dt":41.34, "num_electronic_substeps":1, "nsteps":25000000, "nprint":2500,
              "prefix":"libra_data", "prefix2":"libra_data",
              "hdf5_output_level":-1, "mem_output_level":3, "txt_output_level":-1,
              "use_compression":0, "compression_level":[0,0,0], "progress_frequency":0.05,
              "ntraj":ntraj, "nstates":nstates}

# General Adiabatic recipe
def load_adiabatic(dyn_general, temp=300.0):
    dyn_general.update({"rep_tdse":1}) # adiabatic representation, wfc
    dyn_general.update({"ham_update_method":1})  # recompute only diabatic Hamiltonian
    dyn_general.update({"ham_transform_method":1}) # diabatic->adiabatic according to internal diagonalization
    dyn_general.update({"rep_force":1} ) # adiabatic
    dyn_general.update({"force_method":1} ) # state-specific  as in the TSH or adiabatic
    dyn_general.update({"time_overlap_method":1}) # explicitly compute it from the wavefunction info
    dyn_general.update({"isNBRA":0}) # no NBRA - Hamiltonians for all trajectories are computed explicitly [ default ]
    dyn_general.update({"Temperature":temp}) #Temperature of the system [ default ]
    dyn_general.update({"electronic_integrator":-1}) # No propagation

# KC-RPMD recipe
def load_kcrpmd(dyn_general, kcrpmd_tst, temp=300.0):
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
    dyn_general.update({"kcrpmd_eta":kcrpmd_tst.eta})
    dyn_general.update({"kcrpmd_gamma":kcrpmd_tst.gammay})
    dyn_general.update({"kcrpmd_gammaKP":0.0})
    dyn_general.update({"kcrpmd_my":kcrpmd_tst.my})
    dyn_general.update({"isNBRA":0}) # no NBRA - Hamiltonians for all trajectories are computed explicitly [ default ]
    dyn_general.update({"Temperature":temp}) #Temperature of the system
    dyn_general.update({"electronic_integrator":-1}) # No propagation

# Langevin thermostat for nuclei
def load_langevin(dyn_general, _ndof):
    dyn_general.update({"ensemble":1}) #NVT
    dyn_general.update({"thermostat_params":{ "thermostat_type":"Langevin"
                                             , "Temperature":dyn_general["Temperature"]
                                             , "nu_therm":0.003 }})
    dyn_general.update({"thermostat_dofs":list(range(_ndof))})

# Load in control parameter recipe
if args.method == 1:
    load_adiabatic(dyn_params, T)
    dyn_params.update({"properties_to_save":["timestep","time","q","p","f","Epot_ave","Ekin_ave","Etot_ave"]})
elif (args.method == 2 or args.method == 3):
    load_kcrpmd(dyn_params, kcrpmd_tst, T)
    dyn_params.update({"properties_to_save":["timestep","time","q","p","f","Epot_ave","Ekin_ave","Etot_ave", "y_aux_var","p_aux_var", "f_aux_var", "ekin_aux_var"]})

# Save control parameters for dynamics as copy before setting thermostat and constraints
_control_params_dynamics = dyn_params.copy()
_control_params_dynamics["dt"] = 0.08268
_control_params_dynamics["nsteps"] = 125000
_control_params_dynamics["nprint"] = 100
with open(pref +  "/_control_params_dynamics.txt", "w") as f:
    f.write(str(_control_params_dynamics))

# Load Langevin thermostat for thermalization
load_langevin(dyn_params, ndof)
if args.fix == 's':
    dyn_params.update({"kcrpmd_gamma":0.003})
    dyn_params.update({"kcrpmd_gammaKP":0.003})
    dyn_params.update({"constrained_dofs":[0]})
    dyn_params.update({"quantum_dofs":[]})
    dyn_params["thermostat_dofs"].pop(0)

# Save control parameters for thermalization as copy
_control_params_thermalization = dyn_params.copy()
with open(pref +  "/_control_params_thermalization.txt", "w") as f:
    f.write(str(_control_params_thermalization))

##########################################################################################
# ======= NOW CREATING AND SAVING INITIAL CONDITIONS FOR THERMALIZATION (PART 2) ======= #
##########################################################################################

# As discussed in the paper, free energies do not depend on the mass of coordinates.
# We choose masses so that the harmoic frequencies match the timestep.
# n_therm is redundancy factor, mass of s coordinate will either follow harmonic V0(s), or gaussian restraint for small beta*K0
n_therm = 1000
if args.sys == 3 and args.hw == -1:
    Kref = K0 * np.exp(-bq * 2.1)
else:
    Kref = K0
ms_therm = [max(ms*ws**2*(n_therm*41.34/(2*np.pi))**2, 2*a*beta*(n_therm*ms*ws**2*(s0-s1)*41.34/(2*np.pi*beta*Kref))**2)]
mj_therm = [mj[i]*wj[i]**2*(n_therm*41.34/(2*np.pi))**2 for i in range(len(wj))]
mq_therm = [mq*wq**2*(n_therm*41.34/(2*np.pi))**2]
mass_therm = ms_therm + mj_therm + mq_therm
force_therm = [4 * mass_therm[i] / (beta**2) for i in range(len(mass_therm))]
my_therm = 500000.0

# Initial conditions and masses of nuclear coordinates for thermalization
_nucl_params = {"q":[sdag] + [0.0]*(ndof-2) + [qhw], "p":[0.0]*ndof, "mass": mass_therm,
                "force_constant":force_therm, "init_type":1, "ntraj":ntraj, "ndof": ndof}

# Initial conditions and mass of auxiliary coordinate for thermalization
_elec_params = {"init_type":0, "nstates":nstates, "rep":1, "istate": 0, "ntraj":ntraj, "ndia":ndia, "nadi":nadi,
                "y_aux_var":[ydag], "p_aux_var":[0.0], "m_aux_var":[my_therm]}

# Saving initial conditions for thermalization
with open(pref +  "/_init_nucl_thermalization.txt", "w") as f:
    f.write(str(_nucl_params))

with open(pref +  "/_init_elec_thermalization.txt", "w") as f:
    f.write(str(_elec_params))

