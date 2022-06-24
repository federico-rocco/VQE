# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 16:20:15 2022

@author: cosmo
"""

import mapping as mp 
class Hamiltonian:
    
    def __init__(self, n_fermions, n_qubits, coeff):    #coeff = 2d array
        self.n_fermions = n_fermions
        self.n_qubits = n_qubits
        self.coeff = coeff

        
    def __call__(self): #h_ij a+_i a_j
        h = 0
        for m, n in [[_m, _n] for _m in range(self.n_qubits) for _n in range(self.n_qubits)]:
            h += mp.one_body(self.n_qubits, m, n, coeff=self.coeff[m][n])                
        h = h.reduce()
        return h

