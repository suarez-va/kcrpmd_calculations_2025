import cmath
import sys
import cmath
import math
import os
import h5py

import numpy as np
import matplotlib.pyplot as plt

if sys.platform=="cygwin":
    from cyglibra_core import *
elif sys.platform=="linux" or sys.platform=="linux2":
    from liblibra_core import *
import util.libutil as comn
from libra_py import units
import libra_py.dynamics.heom.compute as compute
from libra_py import ft


from scipy.linalg import inv, pinv, lstsq, svd, eig, hankel
from scipy.integrate import trapezoid, quad
from scipy.special import iv # modified Bessel function of the first kind
from scipy.fft import fft, ifft, fftfreq, fftshift

#os.environ["OMP_NUM_THREADS"] = "8"

plt.rcParams.update({
    'figure.figsize': (6.0, 4.0),
    'figure.dpi': 150,
    'figure.facecolor': 'white',
    'figure.edgecolor': 'white',
    'lines.linewidth': 2,
    'axes.linewidth': 3,
    'axes.labelsize': 15,
    'axes.titlesize': 15,
    'xtick.direction': 'in',
    'xtick.top': True,
    'ytick.direction': 'in',
    'ytick.right': True,
    'xtick.major.width': 1.5,
    'ytick.major.width': 1.5,
    'xtick.major.size': 4,
    'ytick.major.size': 4,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'legend.fontsize': 11,
    'legend.frameon': False,
})

##############################
### Hamiltonian Parameters ###
##############################
#T = 315774.6670514614
T = 300.0
beta = units.hartree / (units.boltzmann * T)
Omega = 0.25 / beta
Lambda = 2.0 / beta
epsilon = 0.0 / beta
Gamma = 0.25 / beta
Delta = 10**(-1.4) / beta


#############################
### Rate Theory Benchmark ###
#############################
P1eq = np.exp(-0.5*beta*epsilon)/(np.exp(-0.5*beta*epsilon)+np.exp(0.5*beta*epsilon))
P2eq = np.exp(0.5*beta*epsilon)/(np.exp(-0.5*beta*epsilon)+np.exp(0.5*beta*epsilon))

kMT12 = 2*np.pi*Delta**2*np.sqrt(beta/(4*np.pi*Lambda))*np.exp(-beta*(Lambda-epsilon)**2/(4*Lambda))
#kNAD12 = 2*np.pi/Omega*Delta**2*np.exp(-0.5*beta*epsilon)*np.exp(-Lambda/(Omega*np.tanh(0.5*beta*Omega)))*iv(-epsilon/Omega, 2*Lambda*np.exp(-0.5*beta*Omega)/(Omega*(1-np.exp(-beta*Omega))))
kNAD12 = kMT12
kZUS12 = kMT12 / (1 + 4*np.pi*Delta**2*Gamma/(Lambda*Omega**2))
#kZUS12 = gam*(beta*Del)**2*(1-(eps/lam)**2)*np.exp(-beta*(lam-eps)**2/(4*lam))/(4*(beta*Del)**2*np.sqrt(np.pi/(beta*lam))+beta*gam*np.sqrt(beta*lam/np.pi)*(1-(eps/lam)**2))

kMT21 = P1eq/P2eq*kMT12
kNAD21 = P1eq/P2eq*kNAD12
kZUS21 = P1eq/P2eq*kZUS12


#####################################
### Spectral Density Distribution ###
#####################################
def wcothm1(w, beta, thr=1e-15):
    is_scalar = np.isscalar(w)
    w = np.atleast_1d(w).astype(float)
    small = np.abs(beta*w) <= thr
    res = np.empty_like(w, dtype=float)
    res[small] = 2/beta - w[small] + (beta*w[small]**2)/6
    res[~small] = -2*w[~small]*np.exp(-beta*w[~small]) / np.expm1(-beta*w[~small])
    if is_scalar:
        return res[0]
    return res

def drude_spectral_density_w(w, Lambda, gamma):
    Jw_w = Lambda*gamma/(gamma**2+w**2)
    return Jw_w

def brownian_spectral_density_w(w, Lambda, Omega, Gamma, Omega_c=None):
    if Omega_c is None:
        Jw_w = Lambda*Omega**2*Gamma/((w**2-Omega**2)**2+Gamma**2*w**2)
    else:
        print("Gotta evaluate Cauchy's Principal Value...")
        Jw_w = 0*w
    return Jw_w

def compute_bath_correlation_function(beta, Jw_w, wpts, dw, tpts, tmax):
    dt = 2*np.pi/(wpts*dw)
    tskip = int(np.floor(tmax/((tpts-1)*dt)))
    tgrid = np.arange(tpts)*tskip*dt
    Ltgrid = np.zeros((tpts), dtype=np.complex128)

    wgrid = np.arange(wpts)*dw
    Jwgrid = Jw_w(wgrid)*wgrid
    Ltgrid = fft(Jwgrid)[np.arange(tpts)*tskip]*dw/np.pi

    Nmin = 0.00001
    Nmax = 750.0
    Ncos = 750.0
    for i, t in enumerate(tgrid):
        func = lambda w: Jw_w(w)*wcothm1(w, beta)*np.cos(w*t)/np.pi
        if t == 0:
            dwint = Nmin/beta
        else:
            dwint = min(Nmin/beta, 2*np.pi/(Ncos*t))
        wintmax = Nmax/beta
        wintpts = int(np.ceil(wintmax/dwint))
        wint = np.linspace(0.0, wintmax, wintpts)
        Ltgrid[i] += trapezoid(func(wint), wint)
        #print(wintpts)
        #wint = np.linspace(0.0, Nwmax/beta, int(np.ceil((Ntmax*t/(2*np.pi))*(Nwmax/beta)))); #print(int(np.ceil((Ntmax*t/(2*np.pi))*(Nwmax/beta))))
        #Ltgrid[i] += quad(func, 0.0, 50/beta, limit=500)[0]
        #if t/beta <= 25:
        #    Ltgrid[i] += quad(func, 0.0, 50/beta, limit=500)[0]
        #    Ltgrid[i] += trapezoid(func(wgrid), wgrid)
        #print(i)
        #func = lambda w: w * Jw_w(w) * np.cos(w*t) / np.pi
        #Ltgrid[i] += trapezoid(func(wgrid), wgrid) # Works
        #func = lambda w: Jw_w(w) * wcothm1(w, beta) * np.cos(w*t) / np.pi
        #Ltgrid[i] += quad(func, 0.0, 25/beta, limit=500)[0]
        #res = quad(func, 0.0, 25/beta, limit=500)[0]
        #print(res)
        #Ltgrid[i] += res
    
    return (wgrid, Jwgrid, tgrid, Ltgrid)

#wpts = 2048000000
wpts = 1024000000
#wpts = 256000000
#wpts = 64000000
#dw = 0.0003125*Omega
#dw = 0.00125*Omega
dw = 0.0025*Omega
#dw = 0.005*Omega
tpts = 1000
tmax = 100.02*beta

Jw_w = lambda w: drude_spectral_density_w(w, Lambda, Omega**2/Gamma)
#Jw_w = lambda w: brownian_spectral_density_w(w, Lambda, Omega, Gamma)

wgrid, Jwgrid, tgrid, Ltgrid = compute_bath_correlation_function(beta, Jw_w, wpts, dw, tpts, tmax)

new = True
#new = False
dt = (tgrid[-1]-tgrid[0])/(tgrid.shape[0]-1)
M=25
print("start")
if new: 
    Ltfit = Ltgrid[1:]
    L=int(np.floor(0.5*(tpts-1)))
    Hankel = hankel(Ltfit[:tpts-1-L], Ltfit[tpts-1-L-1:])
    U, S, W = svd(Hankel); #print(S[:M]/S[0])
    WM0 = W[:M,:L]
    WM1 = W[:M,1:L+1]
    AM = pinv(WM0.T) @ WM1.T
    zM = eig(AM, left=False, right=False)
    VNM = np.vander(zM, tpts-1, increasing=True).T
    aM = -np.log(zM) / dt
    cM = lstsq(VNM, Ltfit)[0]*np.exp(aM*dt)
else:
    L=int(np.floor(0.5*tpts))
    Hankel = hankel(Ltgrid[:tpts-L], Ltgrid[tpts-L-1:])
    U, S, W = svd(Hankel); #print(S[:M]/S[0])
    WM0 = W[:M,:L]
    WM1 = W[:M,1:L+1]
    AM = pinv(WM0.T) @ WM1.T
    zM = eig(AM, left=False, right=False)
    VNM = np.vander(zM, tpts, increasing=True).T
    aM = -np.log(zM) / dt
    cM = lstsq(VNM, Ltgrid)[0]
print("stop")
idx = np.argsort(aM.real)
aM = aM[idx]; cM = cM[idx]

gamma_matsubara = complexList(); c_matsubara = complexList()
setup_bath(M, Lambda, Omega**2/Gamma, T, gamma_matsubara, c_matsubara)
Ltesprit = np.zeros((tpts), dtype=np.complex128)
Ltdrude = np.zeros((tpts), dtype=np.complex128)
for m in range(M):
    Ltesprit += cM[m] * np.exp(-aM[m]*tgrid)
    Ltdrude += c_matsubara[m] * np.exp(-gamma_matsubara[m]*tgrid)
#Error = trapezoid(np.abs((Ltgrid - Ltesprit)/Ltesprit[0]), tgrid) / (tgrid[-1] - tgrid[0])
Error = trapezoid(np.abs((Ltdrude - Ltesprit)/Ltesprit[0]), tgrid) / (tgrid[-1] - tgrid[0])
print(f"δL = {Error}")

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(6, 4.5))
ax1.set_xlabel('$\\beta\hbar\omega$')
ax1.set_ylabel('$J(\omega)/ \Lambda$')
ax1.set_xlim([0.0, 10*beta*Omega**2/Gamma])
ax1.axhline(y = 0.0, linestyle = '--', color = 'k')
ax1.plot(beta*wgrid, Jwgrid/Lambda, color = 'g')

ax2.set_xlabel('$t/\\beta\hbar$')
ax2.set_ylabel('$\\beta\hbar L(t)/ \Lambda$')
ax2.axhline(y = 0.0, linestyle = '--', color = 'k')
#ax2.plot(tgrid/beta, beta*Ltgrid.real/Lambda, color = 'r', label = 'real')
#ax2.plot(tgrid/beta, beta*Ltgrid.imag/Lambda, color = 'b', label = 'imag')

#ax2.plot(tgrid/beta, beta*(Ltgrid - Ltesprit).real/Lambda, color = 'r', label = 'real')
#ax2.plot(tgrid/beta, beta*(Ltgrid - Ltesprit).imag/Lambda, color = 'b', label = 'imag')

#ax2.plot(tgrid/beta, beta*Ltesprit.real/Lambda, color = 'k', linestyle = '--')
#ax2.plot(tgrid/beta, beta*Ltesprit.imag/Lambda, color = 'k', linestyle = '--')

#ax2.plot(tgrid/beta, beta*Ltdrude.real/Lambda, color = 'k', linestyle = '--')
#ax2.plot(tgrid/beta, beta*Ltdrude.imag/Lambda, color = 'k', linestyle = '--')

#ax2.plot(tgrid/beta, beta*(Ltdrude - Ltesprit).real/Lambda, color = 'r', label = 'real')
#ax2.plot(tgrid/beta, beta*(Ltdrude - Ltesprit).imag/Lambda, color = 'b', label = 'imag')
ax2.set_ylim([-1e-7,1e-7])

ax2.plot(tgrid/beta, beta*(Ltgrid - Ltdrude).real/Lambda, color = 'r', label = 'real')
ax2.plot(tgrid/beta, beta*(Ltgrid - Ltdrude).imag/Lambda, color = 'g', label = 'imag')
#ax2.set_ylim([-1e-8,1e-8])
#ax2.plot(tgrid/beta, beta*(Ltgrid.imag + 0.5*Lambda*Omega**2/Gamma*np.exp(-Omega**2/Gamma*tgrid))/Lambda, color = 'b', label = 'imag')
#print(gamma_matsubara[0], Omega**2/Gamma)
#print(c_matsubara[0], 0.5*Lambda*Omega**2/Gamma*(1/np.tan(0.5*beta*Omega**2/Gamma) - 1j))

ax2.legend(ncol=2, fontsize=9, loc='best', handletextpad=0.2,  columnspacing=0.9)
plt.subplots_adjust(hspace=0.4)
plt.show()
#plt.close()

#gamma_matsubara = complexList(); c_matsubara = complexList()
#setup_bath(M, Lambda, Omega**2/Gamma, T, gamma_matsubara, c_matsubara)
for m in range(M):
    print(f"gamma_matsubara.append{aM[m]}; c_matsubara.append{cM[m]}")
    #print(f"gamma_matsubara.append{aM[-(1+m)]}; c_matsubara.append{cM[-(1+m)]}")
    print(f"gamma_matsubara.append{gamma_matsubara[m]}; c_matsubara.append{c_matsubara[m]}")

print(f"{Jwgrid.nbytes / 1024**2:.2f} MB")

