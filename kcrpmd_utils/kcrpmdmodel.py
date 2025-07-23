import os
import sys
import math
import copy

if sys.platform == "cygwin":
    from cyglibra_core import *
elif sys.platform == "linux" or sys.platform == "linux2":
    from liblibra_core import *
import util.libutil as comn
import libra_py.units as units

import numpy as np
from scipy.optimize import minimize

class tmp:
    pass

def gen_kcrpmd_bath_params(_params):
    """
    Generates the parameters

    Args:
        params ( dictionary ): model parameters, can contain:

            * **params["M"]** ( double ): Ohmic bath coordinate mass(es) [ default: 1836., units: a.u. ]
            * **params["wc"]** ( double ): Ohmic bath frequency [ default: 2.28e-3, units: a.u. ]
            * **params["gam"]** ( double ): Ohmic bath strength [ default: 4.18608, units: a.u. ]
            * **params["f"]** ( integer ): number of Ohmic bath dof [ default: 0 ]

    Returns:
        tuple: (list, list, list)

            * omega ( list of `num_osc` doubles ): bath frequencies [ units: Ha/Bohr^2 ]
            * coupl ( list of `num_osc` doubles ): coupling strengths [ units: Ha/Bohr ]
            * mass ( list of `num_osc` doubles): masses of oscillators [ units: a.u. ]
    """

    params = dict(_params)

    critical_params = []
    default_params = {
        "M": 1836,
        "wc": 2.28e-3,
        "gam": 4.18608,
        "f": 0}
    comn.check_input(params, default_params, critical_params)

    M = params["M"]
    wc = params["wc"]
    gam = params["gam"]
    f = params["f"]

    omega, coupl, mass = [], [], []

    # Ohmic bath discretizations are according to:
    # https://doi.org/10.1063/1.4863919

    pref = math.sqrt(2 * gam * M * wc / (f * math.pi))
    omega = [-wc * math.log((k - 0.5) / f) for k in range(1, f + 1)]
    coupl = [omega[k] * pref for k in range(f)]
    mass = [M for k in range(f)]

    return omega, coupl, mass


def get_ABC(q0, leps, Ea):
    def minimize_me(x):
        A, B, C = x
        A_abs, B_abs = np.abs(A), np.abs(B)

        discriminant = 9 * B**2 - 32 * A_abs * C
        if discriminant < 0:
            return 1e12  # large penalty to avoid invalid region

        delta = np.sqrt(discriminant)
        q_plus = (3 * B_abs + delta) / (8 * A_abs)
        q_minus = (3 * B_abs - delta) / (8 * A_abs)

        E_plus = A_abs * q_plus**4 - B_abs * q_plus**3 + C * q_plus**2
        E_minus = A_abs * q_minus**4 - B_abs * q_minus**3 + C * q_minus**2

        term_q0 = (q0 - q_plus)**2
        term_leps = (leps - E_plus)**2
        term_Ea = (Ea - E_minus)**2

        return term_q0 + term_leps + term_Ea

    x_guess = np.array([16 * Ea / q0**4, 32 * Ea / q0**3, 16 * Ea / q0**2])

    for _ in range(7):
        res = minimize(minimize_me, x_guess, method='Nelder-Mead', tol=1e-11)
        x_guess = res.x

    Aq, Bq, Cq = float(np.abs(x_guess[0])), float(-np.abs(x_guess[1])), float(x_guess[2])
    return Aq, Bq, Cq


def kcrpmd_system_bath(q, params, full_id):
    """
   
    effective solvent coordinate factored 2-state spin-boson model as defined in KC-RPMD paper

         | 0.5*ms*ws^2*(s-s0)^2            K0*exp(-bq*q)       |
    H =  |                                                     | + I*[Vb(s,x) + Vda(q)]
         |    K0*exp(-bq*q)         0.5*ms*ws^2*(s-s1)^2 + eps |
         

    Vb(s,x) = sumj{0.5*Mj*wj^2*(xj-cj*s/(Mj*wj^2))^2}

    Ohmic spectral density as defined in KC-RPMD paper (gamma = pi *hbar / 2 * ksi, ksi - Kondo parameter)
    J(w) = gam * w * exp(-|w|/wc)

    Args:
        q ( MATRIX(ndof, 1) ): coordinates of the particle, ndof = f
        params ( dictionary ): model parameters
            * **params["ms"]** ( double ): s coordinate mass [ default: 1836., units: a.u. ]
            * **params["ws"]** ( double ): s coordinate angular frequency [ default: 2.28e-3, units: a.u. ]
            * **params["s0"]** ( double ): V0 parabola center [ default: -2.40, units: a.u. ]
            * **params["s1"]** ( double ): V1 parabola center [ default: 2.40, units: a.u. ]
            * **params["eps"]** ( double ): V1 parabola vertical shift [ default: 0., units: a.u. ]
            * **params["wj"]** (list of f doubles): frequencies [ default: [], units: a.u.]
            * **params["cj"]** (list of f doubles): diagonal linear couplings [ default: [], units: a.u.]
            * **params["Mj"]** (list of f doubles): masses of nuclear DOFS [ default: [], units: a.u.]
            * **params["K0"]** ( double ): electronic coupling strength [ default: 6.67e-7, units: a.u. ]
            * **params["sys_type"]** ( int ): Different options for Vda(q):
              - 0: do not include the donor acceptor coordinate q [ default ]
              - 1: system A of KC-RPMD paper
              - 2: system B of KC-RPMD paper
              - 3: system C of KC-RPMD paper
            * **params["bq"]** ( double ): exponential coupling parameter [ default: 0., units: a.u. ]
            * **params["mq"]** ( double ): q coordinate mass [ default: 5.00e4, units: a.u. ]
            * **params["wq"]** ( double ): systems A and B q coordinate angular frequency [ default: 5.00e-4, units: a.u. ]
            * **params["Dq"]** ( double ): systems B and C morse potential dissociation [ default: 1.00e-3, units: a.u. ]
            * **params["Aq"]** ( double ): system C quartic coefficient [ default: 1.041e-2, units: a.u. ]
            * **params["Bq"]** ( double ): system C cubic coefficient [ default: 4.065e-2, units: a.u. ]
            * **params["Cq"]** ( double ): system C quadratic coefficient [ default: 3.622e-2, units: a.u. ]
            * **params["hard_wall"]** ( int ): whether to set a sextic hard wall potential for Vda(q):
              - -1: left side hard wall
              - 0: no hard wall [ default ]
              - 1: right side hard wall
            * **params["qhw"]** ( double ): hard wall position [ default: 1.00, units: a.u. ]
            * **params["khw"]** ( double ): hard wall constant [ default: 1.00e5, units: a.u. ]

    Returns:       
        PyObject: obj, with the members:

            * obj.ham_dia ( CMATRIX(2,2) ): diabatic Hamiltonian
            * obj.ovlp_dia ( CMATRIX(2,2) ): overlap of the basis (diabatic) states [ identity ]
            * obj.d1ham_dia ( list of ndof CMATRIX(2,2) objects ):
                derivatives of the diabatic Hamiltonian w.r.t. the nuclear coordinate
            * obj.dc1_dia ( list of ndof CMATRIX(2,2) objects ): derivative coupling in the diabatic basis [ zero ]
 
    """

    critical_params = [] 
    default_params = {"ms":1836, "ws":2.28e-3, "s0":-2.4, "s1":2.4, "eps":0,
                      "wj":[], "cj":[], "Mj":[], "K0":6.67e-7,
                      "sys_type":0, "bq":0., "mq":5e4, "wq":5e-4,
                      "Dq":1e-3, "Aq":1.041e-2, "Bq":-4.065e-2, "Cq":3.622e-2,
                      "hard_wall":0, "qhw":1., "khw":1e5}
    comn.check_input(params, default_params, critical_params)

    ms = params["ms"]
    ws = params["ws"]
    s0 = params["s0"]
    s1 = params["s1"]
    eps = params["eps"]
    wj = params["wj"]
    cj = params["cj"]
    Mj = params["Mj"]
    K0 = params["K0"]
    sys_type = params["sys_type"]
    bq = params["bq"]
    mq = params["mq"]
    wq = params["wq"]
    Dq = params["Dq"]
    Aq = params["Aq"]
    Bq = params["Bq"]
    Cq = params["Cq"]
    hard_wall = params["hard_wall"]
    qhw = params["qhw"]
    khw = params["khw"]

    f = len(Mj)
    ndof = 2 + f if (sys_type != 0) else 1 + f
    if q.num_of_rows != ndof:
        print(f"Shape of coordinates inconsistent with system parameters, q should have {ndof} rows\nExiting now...\n")
        sys.exit(0)

    obj = tmp()
    obj.ham_dia = CMATRIX(2, 2)
    obj.ovlp_dia = CMATRIX(2, 2);  obj.ovlp_dia.identity()
    obj.d1ham_dia = CMATRIXList()
    obj.dc1_dia = CMATRIXList()

    for i in range(ndof):
        obj.d1ham_dia.append( CMATRIX(2, 2) )
        obj.dc1_dia.append( CMATRIX(2, 2) )

    indx = 0
    if full_id !=None:
        Id = Cpp2Py(full_id)
        indx = Id[-1]

    #=========== Energies & Derivatives ===============

    # s is the effective solvent coordinate, unique from other bath coordinates x
    s = q.get(0, indx)

    # energy from just s coordinate without bath x coordinates
    obj.ham_dia.set(0, 0, 0.5 * ms * ws**2 * (s - s0)**2 * (1.0+0.0j))
    obj.ham_dia.set(1, 1, (0.5 * ms * ws**2 * (s - s1)**2 - eps) * (1.0+0.0j))

    # derivative w.r.t. s from just s coordinate without bath x coordinates
    obj.d1ham_dia[0].add(0, 0, ms * ws**2 * (s - s0) * (1.0+0.0j))
    obj.d1ham_dia[0].add(1, 1, ms * ws**2 * (s - s1) * (1.0+0.0j))

    x = 0.0
    for j in range(1, f + 1):
        x_j = q.get(j, indx)

        # energy

        c_j = cj[j-1]
        w_j = Mj[j-1] * wj[j-1]**2
        u_j = x_j - c_j * s / w_j

        x += 0.5 * w_j * u_j**2
        y = -c_j * u_j
        z = w_j * u_j

        # derivative w.r.t. s:
        obj.d1ham_dia[0].add(0, 0, y * (1.0+0.0j))
        obj.d1ham_dia[0].add(1, 1, y * (1.0+0.0j))
        # derivative w.r.t. x_j:
        obj.d1ham_dia[j].add(0, 0, z * (1.0+0.0j))
        obj.d1ham_dia[j].add(1, 1, z * (1.0+0.0j))
        
    obj.ham_dia.add(0, 0, x * (1.0+0.0j))
    obj.ham_dia.add(1, 1, x * (1.0+0.0j))

    # qda is the donor-acceptor coordinate
    if (sys_type == 0):
        obj.ham_dia.set(0, 1, K0 * (1.0+0.0j))
        obj.ham_dia.set(1, 0, K0 * (1.0+0.0j))
    else:
        qda = q.get(f + 1, indx)
        obj.ham_dia.set(0, 1, K0 * np.exp(-bq * qda) * (1.0+0.0j))
        obj.ham_dia.set(1, 0, K0 * np.exp(-bq * qda) * (1.0+0.0j))
        obj.d1ham_dia[f + 1].set(0, 1, -bq * K0 * np.exp(-bq * qda) * (1.0+0.0j))
        obj.d1ham_dia[f + 1].set(1, 0, -bq * K0 * np.exp(-bq * qda) * (1.0+0.0j))
        if ((hard_wall == -1 and qda <= qhw) or (hard_wall == 1 and qda >= qhw)):
            obj.ham_dia.add(0, 0, khw * (qda - qhw)**6 * (1.0+0.0j))
            obj.ham_dia.add(1, 1, khw * (qda - qhw)**6 * (1.0+0.0j))
            obj.d1ham_dia[f + 1].add(0, 0, 6 * khw * (qda - qhw)**5 * (1.0+0.0j))
            obj.d1ham_dia[f + 1].add(1, 1, 6 * khw * (qda - qhw)**5 * (1.0+0.0j))
        if (sys_type == 1):
            obj.ham_dia.add(0, 0, 0.5 * mq * wq**2 * qda**2 * (1.0+0.0j))
            obj.ham_dia.add(1, 1, 0.5 * mq * wq**2 * qda**2 * (1.0+0.0j))
            obj.d1ham_dia[f + 1].add(0, 0, mq * wq**2 * qda * (1.0+0.0j))
            obj.d1ham_dia[f + 1].add(1, 1, mq * wq**2 * qda * (1.0+0.0j))
        elif (sys_type == 2):
            if (qda > 0.0):
                obj.ham_dia.add(0, 0, 0.5 * mq * wq**2 * qda**2 * (1.0+0.0j))
                obj.ham_dia.add(1, 1, 0.5 * mq * wq**2 * qda**2 * (1.0+0.0j))
                obj.d1ham_dia[f + 1].add(0, 0, mq * wq**2 * qda * (1.0+0.0j))
                obj.d1ham_dia[f + 1].add(1, 1, mq * wq**2 * qda * (1.0+0.0j))
            else:
                obj.ham_dia.add(0, 0, Dq * (1 - np.exp(-np.sqrt(0.5 * mq * wq**2 / Dq) * qda))**2 * (1.0+0.0j))
                obj.ham_dia.add(1, 1, Dq * (1 - np.exp(-np.sqrt(0.5 * mq * wq**2 / Dq) * qda))**2 * (1.0+0.0j))
                obj.d1ham_dia[f + 1].add(0, 0, 2 * Dq * np.sqrt(0.5 * mq * wq**2 / Dq) * np.exp(-np.sqrt(0.5 * mq * wq**2 / Dq) * qda) * (1 - np.exp(-np.sqrt(0.5 * mq * wq**2 / Dq) * qda)) * (1.0+0.0j))
                obj.d1ham_dia[f + 1].add(1, 1, 2 * Dq * np.sqrt(0.5 * mq * wq**2 / Dq) * np.exp(-np.sqrt(0.5 * mq * wq**2 / Dq) * qda) * (1 - np.exp(-np.sqrt(0.5 * mq * wq**2 / Dq) * qda)) * (1.0+0.0j))
        elif (sys_type == 3):
            if (qda > 0.0):
                obj.ham_dia.add(0, 0, Aq * qda**4 + Bq * qda**3 + Cq * qda**2 * (1.0+0.0j))
                obj.ham_dia.add(1, 1, Aq * qda**4 + Bq * qda**3 + Cq * qda**2 * (1.0+0.0j))
                obj.d1ham_dia[f + 1].add(0, 0, 4 * Aq * qda**3 + 3 * Bq * qda**2 + 2 * Cq * qda * (1.0+0.0j))
                obj.d1ham_dia[f + 1].add(1, 1, 4 * Aq * qda**3 + 3 * Bq * qda**2 + 2 * Cq * qda * (1.0+0.0j))
            else:
                obj.ham_dia.add(0, 0, Dq * (1 - np.exp(-np.sqrt(Cq / Dq) * qda))**2 * (1.0+0.0j))
                obj.ham_dia.add(1, 1, Dq * (1 - np.exp(-np.sqrt(Cq / Dq) * qda))**2 * (1.0+0.0j))
                obj.d1ham_dia[f + 1].add(0, 0, 2 * Dq * np.sqrt(Cq / Dq) * np.exp(-np.sqrt(Cq / Dq) * qda) * (1 - np.exp(-np.sqrt(Cq / Dq) * qda)) * (1.0+0.0j))
                obj.d1ham_dia[f + 1].add(1, 1, 2 * Dq * np.sqrt(Cq / Dq) * np.exp(-np.sqrt(Cq / Dq) * qda) * (1 - np.exp(-np.sqrt(Cq / Dq) * qda)) * (1.0+0.0j))
    
    return obj

