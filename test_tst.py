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

# ======= TST code to evaluate eta, gamma, and mass of auxiliary variable =======
T = 300.0 
beta = units.hartree / (units.boltzmann * T)
a = 0.1
b = 1000.0
c = 0.5
d = 3.0

ms = 1836.0
ws = 2.28e-3
s0 = -2.4 
s1 = 2.4 
eps = 0.0

K0 = 2.85e-3
bq = 3.0
Dq = 1.0e-4
(q0, Ea, leps) = (2.1, 6.65e-3, -1.62e-2)
(Aq, Bq, Cq) = get_ABC(q0, leps, Ea)
qhw = 1.0
khw = 1.0e5

Kq = lambda q: K0 * np.exp(-bq * q)
Vq = lambda q: np.piecewise(q, [q >= 0., q < 0.],
                            [lambda q: Aq * q**4 + Bq * q**3 + Cq * q**2,
                             lambda q: Dq * (1 - np.exp(-np.sqrt(Cq / Dq) * q))**2])

kcrpmd_tst = KcrpmdTst(beta, a, b, c, d, 1., ms, ws, s0, s1, eps, Kq, Vq)
kcrpmd_tst.set_eta_my_gammay()
print(kcrpmd_tst.eta, kcrpmd_tst.my)
#kGR = kcrpmd_tst.kGR()
#ktsty = kcrpmd_tst.tst_y()
#print(kGR, ktsty)

print("progress 1")
hw = -1
if hw != 0:
    q_low_cp = kcrpmd_tst.q_low
    print(q_low_cp)
    kcrpmd_tst.q_low = qhw
    print(kcrpmd_tst.q_low)
    kcrpmd_tst.set_eta_my_gammay()
    kcrpmd_tst.q_low = q_low_cp
    print(q_low_cp)
    print(kcrpmd_tst.q_low)

print(kcrpmd_tst.eta, kcrpmd_tst.my)
kGR = kcrpmd_tst.kGR()
ktsty = kcrpmd_tst.tst_y()
print(kGR, ktsty)

exit()


kcrpmd_tst.set_eta_my_gammay()
print(kcrpmd_tst.eta, kcrpmd_tst.my)
kcrpmd_tst.q_low = qhw
kcrpmd_tst.set_eta_my_gammay()
print(kcrpmd_tst.eta, kcrpmd_tst.my)
kGR = kcrpmd_tst.kGR()
ktsty = kcrpmd_tst.tst_y()
print(kGR, ktsty)


if hw == -1:
    kcrpmd_tst.q_low = qhw - (100 / (khw * beta))**(1/6)
    Vhw = lambda q: np.piecewise(q, [q >= qhw, q < qhw], [lambda q: np.full_like(q, 0.), lambda q: khw * (q - qhw)**6])
    kcrpmd_tst.Vq = lambda q: Vq(q) + Vhw(q) 
elif hw == 1:
    kcrpmd_tst.q_high = qhw + (100 / (khw * beta))**(1/6)
    Vhw = lambda q: np.piecewise(q, [q >= qhw, q < qhw], [lambda q: khw * (q - qhw)**6, lambda q: np.full_like(q, 0.)])
    kcrpmd_tst.Vq = lambda q: Vq(q) + Vhw(q) 

print(kcrpmd_tst.eta, kcrpmd_tst.my)
kGR = kcrpmd_tst.kGR()
ktsty = kcrpmd_tst.tst_y()
print(kGR, ktsty)
kcrpmd_tst.set_eta_my_gammay()
print(kcrpmd_tst.eta, kcrpmd_tst.my)
kGR = kcrpmd_tst.kGR()
ktsty = kcrpmd_tst.tst_y()
print(kGR, ktsty)

