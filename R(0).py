# -*- coding: utf-8 -*-
"""
Created on Mon Oct 17 10:53:31 2022

@author: cosmo
"""

import numpy as np
import math as mt
import scipy as sp
b = 0.0016693380059727566*1000
def psi(n):
    return (-1)**n*np.sqrt(2*mt.factorial(n)/(b**3*mt.gamma(n+3/2)))*sp.special.genlaguerre(n, 1/2)(0)

print((-0.9858*psi(1)-0.1369*psi(2)-0.09722*psi(3))**2)