# -*- coding: utf-8 -*-
"""
Created on Sat Jun 11 16:02:16 2022

@author: cosmo
"""
from ferm_op import Hamiltonian
from qiskit.aqua.algorithms import NumPyMinimumEigensolver
import numpy as np
import scipy as sp
import mpmath
from qiskit.opflow import I, X, Z, Y

k=0.4063
sigma=441.6**2
mu=637.5
omega=562.9
b=1/np.sqrt(mu*omega)

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

# Instantiate the system Hamiltonian
hamiltonian = (1/2*(21/4*omega + V(0,0) + V(1,1) + V(2,2)) * I ^ I ^ I) + \
        (-1/2*(3/4*omega + V(0,0)) * I ^ I ^ Z) + \
        (-1/2*(7/4*omega + V(1,1)) * I ^ Z ^ I) + \
        (-1/2*(11/4*omega + V(2,2)) * Z ^ I ^ I) + \
        (1/4*(-np.sqrt(3/2)*omega + 2*V(0,1)) * I ^ X ^ X) + \
        (1/4*(-np.sqrt(5)*omega + 2*V(1,2)) * X ^ X ^ I) + \
        (1/4*(-np.sqrt(3/2)*omega + 2*V(0,1)) * I ^ Y ^ Y) + \
        (1/4*(-np.sqrt(5)*omega + 2*V(1,2)) * Y ^ Y ^ I) + \
        (1/2*V(0,2) * X ^ Z ^ X) + \
        (1/2*V(0,2) * Y ^ Z ^ Y) 
nm = NumPyMinimumEigensolver(operator=hamiltonian)
print(nm.compute_minimum_eigenvalue())