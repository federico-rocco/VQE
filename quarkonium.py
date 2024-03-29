# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 18:47:14 2022

@author: cosmo
"""
import numpy as np
import scipy as sp
import mpmath
quarkonium_parameters = { #mass, #mass for spin, number of flavours, sigma, oscillation frequency, |R(0)|**2
    'charmonium'   :[ 1317, 1317, 1410, 4,  0.18*10**6, 562.9, 0.856],
    'bottomonium'  :[ 4584,  4584, 4769, 5,  0.25*10**6, 400, 6.476],
    'bc'           :[ 1317,  4584, 4769, 5,  0.25*10**6, 500, 6.476]
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
        self.mq, self.mqbar, self.mqs, nf, self.sigma, self.omega, self.R0 = quarkonium_parameters[flavour]
        μ = 2*self.mq*self.mqbar/(self.mq+self.mqbar) #energy scale
        self.αs = 12*np.pi/((33-2*nf)*np.log(μ**2/Λ**2))
        self.k = (4/3)*self.αs
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
        return 1000*8*αHF*self.R0/(9*(self.mqs/1000)**2)*(S*(S+1)/2-3/4) 
    
    
    def coeff(self, orbitals):
        coeff = np.zeros((orbitals, orbitals))
        for m, n in [[_m, _n] for _m in range(orbitals) for _n in range(orbitals)]:
            coeff[m][n] += self.T(m,n) - self.k*self.r_inv(m,n) + self.sigma*self.r(m,n)
        return [coeff,None]
    
