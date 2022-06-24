# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 16:20:15 2022

@author: cosmo
"""

import mapping as mp 
import numpy as np
class Hamiltonian:
    
    def __init__(self, n_fermions, n_qubits, coeff):    #coeff = 2d array
        self.n_fermions = n_fermions
        self.n_qubits = n_qubits
        self.coeff = coeff
        self.list = []

        
    def __call__(self): #h_ij a+_i a_j
        h = 0
        for m, n in [[_m, _n] for _m in range(self.n_qubits) for _n in range(self.n_qubits)]:
             h +=mp.one_body(self.n_qubits, m, n, coeff=self.coeff[m][n])            
        h = h.reduce()
        
        pauli_list = []
        for pauli_sum_op in h:
            pauli_op = pauli_sum_op.to_pauli_op()
            pauli_list.append(Pauli(pauli_op))

        return pauli_list

class Pauli: #PauliOp with metods
    
    def __init__(self, pauli):
        self.pauli = pauli #PauliOp object
        self.pauli_string = str(pauli.primitive) #corresponding string
        self.coeff = pauli.coeff
        
    def pauli_to_qc(self, qc, qb):
        #print(self.pauli_string)
        for i, item in enumerate(self.pauli_string):
            #print(i,item)
            if item == 'X':
                qc.h(qb[i])
                #self.pauli_string[i] = 'Z'
            elif item == 'Y':
                qc.rx(-np.pi/2, qb[i])
                #self.pauli_string[i] = 'Z'
        #print("string applied, " ,self.pauli_string)
        return qc
    
    def expectation(self, measurement, shots):
        #get expectation value
        
        exp_value = 0
        #print(self.pauli_string)
        for (state,value) in measurement.items():
            #print(state)
            sign = -1
            for (i,number) in enumerate(state):
                if number == 1 and self.pauli_string[i] != 'I':
                    #print("Z and 1 in position ",i)
                    sign *= -1 
            exp_value += sign * value/shots
        return exp_value
    
    #def z_count(self):
        #for 
    
    
    
    
    
    
    
    
    
    
    
    