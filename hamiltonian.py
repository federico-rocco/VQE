# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 16:20:15 2022

@author: cosmo
"""

import mapping as mp 
import numpy as np

import qiskit


class Hamiltonian:
    
    def __init__(self, n_fermions=None, n_qubits=None, coeff=None):    #coeff = 2d array
        self.n_fermions = n_fermions
        self.n_qubits = n_qubits
        self.coeff = coeff

        
    def __call__(self, ising=False): #h_ij*a+_i*a_j
    
        if ising == True:
            return self.isingH()
        
        h = 0
        for m, n in [[_m, _n] for _m in range(self.n_qubits) for _n in range(self.n_qubits)]:
             h +=mp.one_body(self.n_qubits, m, n, coeff=self.coeff[m][n])            
        h = h.reduce()
        
        pauli_list = []
        for pauli_sum_op in h:
            pauli_op = pauli_sum_op.to_pauli_op()
            pauli_list.append(Pauli(pauli_op))

        
        return pauli_list
    
    def isingH(self):
        pauli_list = []
        
        #nearest neighbours term
        for i in range(self.n_qubits):
            label_nn = ""
            for item in range(self.n_qubits):                 
                if item == i:
                    label_nn += "Z"
                else:
                    label_nn += "I"
                    
            op_nn = Pauli(qiskit.opflow.PauliOp(qiskit.quantum_info.Pauli(label_nn)))
            op_nn.coeff = -self.coeff['B']/2
            pauli_list.append(op_nn)
            
         
        #lattice term
        for i in range(self.n_qubits):
            label_lat = ""
            
            if i == self.n_qubits-1:
                print("hey")
                for item in range(self.n_qubits):                 
                    if item == i or item == 0:
                        label_lat += "Z"
                    else:
                        label_lat += "I"
            
            else:
                for item in range(self.n_qubits):                 
                    if item == i or item == i+1:
                        label_lat += "Z"
                    else:
                        label_lat += "I"
                    
            op_lat = Pauli(qiskit.opflow.PauliOp(qiskit.quantum_info.Pauli(label_lat)))
            op_lat.coeff = -self.coeff['J']/4
            pauli_list.append(op_lat)
            
        
        for p in pauli_list:
            print(p.pauli_string)
        return pauli_list
            

class Pauli: #PauliOp with metods
    
    def __init__(self, pauli):
        self.pauli = pauli #PauliOp object
        self.pauli_string = str(pauli.primitive) #corresponding string
        self.coeff = pauli.coeff
        
    def pauli_to_qc(self, qc, qb):
        #basis change 
        for i, item in enumerate(self.pauli_string):
            if item == 'X':
                qc.h(qb[i])
            elif item == 'Y':
                qc.rx(-np.pi/2, qb[i])
        return qc
    
    def expectation(self, measurement, shots):
        #get expectation value
        #print(measurement)
        exp_value = 0
        for (state,value) in measurement.items():
            sign = 1
            for (i,number) in enumerate(state):
                if (number == '1' and self.pauli_string[i] != 'I'):
                    sign *= -1 
            exp_value += sign * value/shots
        return exp_value

    
    
    
    
    
    
    
    
    
    
    
    