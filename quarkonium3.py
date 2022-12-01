# -*- coding: utf-8 -*-
"""
Created on Thu Oct 20 11:19:54 2022

@author: cosmo
"""

import numpy as np
import scipy as sp
import mpmath
quarkonium_parameters = { #mass, #mass for spin, number of flavours, sigma, oscillation frequency, |R(0)|**2
    'charmonium'   :[ 1240, 0.4,  0.2355*10**6, 562.9, 0.662],
    'bottomonium'  :[ 4500,  0.3,  0.4237*10**6, 400, 6.170],
    'bc'           :[ (1240+4500)/2,  0.34,  0.2986*10**6, 500, 1.381]
}
Λ = 150 #Λqcd
αHF = 0.31

def kron_delta(a,b):
    if a == b :
        return 1
    else:
        return 0
    
class Quarkonium:
    
    def __init__(self, flavour):
        self.mq, self.αs, self.sigma, self.omega, self.R0 = quarkonium_parameters[flavour]
        μ = self.mq #energy scale
        self.k = self.αs
        mu = μ/2 #reduced mass
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
    
    
    def VSS(self, S):
        return 1000*8*self.αs*self.R0/(9*(self.mq/1000)**2)*(S*(S+1)/2-3/4) 
    
    
    def coeff(self, orbitals):
        coeff = np.zeros((orbitals, orbitals))
        for m, n in [[_m, _n] for _m in range(orbitals) for _n in range(orbitals)]:
            coeff[m][n] += self.T(m,n) - self.k*self.r_inv(m,n) + self.sigma*self.r(m,n)
        return [coeff,None]