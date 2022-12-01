# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 09:59:25 2022

@author: cosmo
"""

import numpy as np
import scipy as sp
import mpmath


def kron_delta(a,b):
    if a == b :
        return 1
    else:
        return 0
    
class Quarkonium:
    
    def __init__(self, flavour='charmonium'):
        self.omega = 562.9
        self.sigma = 441.6**2
        self.k = 0.4063
        mu = 637.5
        self.b = 1/np.sqrt(mu*self.omega)   
    
    
    def T(self, m, n):
        return (self.omega/2*((2*n+3/2)*kron_delta(m,n)-
                   np.sqrt(n*(n+1/2))*kron_delta(m+1,n)-
                   np.sqrt((n+1)*(n+3/2))*kron_delta(m-1,n)))
        
    def r(self, m, n):
        return ((-1)**(m+n)*4*self.b/(np.pi*(1-4*n**2))*np.sqrt(sp.special.gamma(m+3/2)*
                   sp.special.gamma(n+3/2)/(sp.special.factorial(m)*sp.special.factorial(n)))*
                   float(mpmath.hyp2f1(2,-m,3/2-n,1)))
        
    def r_inv(self, m, n):
        return ((-1)**(m+n)*4/(self.b*np.pi*(1+2*n))*np.sqrt(sp.special.gamma(m+3/2)*
                   sp.special.gamma(n+3/2)/(sp.special.factorial(m)*sp.special.factorial(n)))*
                   float(mpmath.hyp3f2(1/2,1,-m,3/2,1/2-n,1)))
    
    
    
    def coeff(self, orbitals):
        coeff = np.zeros((orbitals, orbitals))
        for m, n in [[_m, _n] for _m in range(orbitals) for _n in range(orbitals)]:
            coeff[m][n] += self.T(m,n) - self.k*self.r_inv(m,n) + self.sigma*self.r(m,n)
        return [coeff,None]