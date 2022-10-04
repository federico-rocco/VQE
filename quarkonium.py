# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 18:47:14 2022

@author: cosmo
"""
import numpy as np
import scipy as sp
import mpmath

#nf = 4 #number of active flavours at μ
μ = 1317 #2(mq*mqbar)/(mq+mqbar)
Λ = 150 #Λqcd
αs = 0.3047#12*np.pi/((33-2*nf)*np.log(μ**2/Λ**2))
#αs2 = 0.220245

#sigma = 441.6**2 #MeV^2#0.18*10**6
mu = 1317/2 #637.5
omega = 562.9
b = 1/np.sqrt(mu*omega)
#mq = 1.4794 #GeV
#mq2 = 1270 #MeV
orbitals = 3
fermions = 1

mc = 1317 #MeV
quarkonium_parameters = { #mass, nf, sigma
    'charmonium'   :[ 1317,  3,  0.18*10**6],
    'bottomonium'  :[ 4584,  4,  0.25*10**6]
}

def kron_delta(a,b):
    if a == b :
        return 1
    else:
        return 0

def T(m,n):
    return (omega/2*((2*n+3/2)*kron_delta(m,n)-
               np.sqrt(n*(n+1/2))*kron_delta(m+1,n)-
               np.sqrt((n+1)*(n+3/2))*kron_delta(m-1,n)))
    
def r(m,n):
    return ((-1)**(m+n)*4*b/(np.pi*(1-4*n**2))*np.sqrt(sp.special.gamma(m+3/2)*
               sp.special.gamma(n+3/2)/(sp.special.factorial(m)*sp.special.factorial(n)))*
               float(mpmath.hyp2f1(2,-m,3/2-n,1)))
    
def r_inv(m,n):
    return ((-1)**(m+n)*4/(b*np.pi*(1+2*n))*np.sqrt(sp.special.gamma(m+3/2)*
               sp.special.gamma(n+3/2)/(sp.special.factorial(m)*sp.special.factorial(n)))*
               float(mpmath.hyp3f2(1/2,1,-m,3/2,1/2-n,1)))


def VSS(S):
    return (32*np.pi*αs/(9*(mc*10**(-3))**2))*(S*(S+1)/2-3/4) 


def coeff(flavour):
    mass, nf, sigma = quarkonium_parameters[flavour]
    alphas = np.pi/((33-2*nf)*np.log(mass**2/Λ**2))
    k = 4/3*0.3047
    coeff = np.zeros((orbitals, orbitals))
    for m, n in [[_m, _n] for _m in range(orbitals) for _n in range(orbitals)]:
        coeff[m][n] += T(m,n) -k*r_inv(m,n) + 441.6**2*r(m,n)
    return [coeff,None]