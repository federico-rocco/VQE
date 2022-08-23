# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 18:47:14 2022

@author: cosmo
"""
import numpy as np
import scipy as sp
import mpmath


k=0.4063 #4/3*αs 
sigma=441.6**2 #MeV^2
mu=637.5
omega=562.9
b=1/np.sqrt(mu*omega)
orbitals = 3
fermions = 1


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
    
def V(m,n):
    return -k*r_inv(m,n)+sigma*r(m,n)


def coeff():
    coeff = np.zeros((orbitals, orbitals))
    for m, n in [[_m, _n] for _m in range(orbitals) for _n in range(orbitals)]:
        coeff[m][n] += T(m,n) + V(m,n)
    return coeff