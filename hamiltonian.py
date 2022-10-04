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

        
    def __call__(self, ising=False): 
    
        if ising == True:
            return self.isingH()
        
        h = 0
        
        #h_ij * a+_i * a_j
        for i, j in [[_i, _j] for _i in range(self.n_qubits) for _j in range(self.n_qubits)]:
             h += mp.one_body(self.n_qubits, i, j, coeff=self.coeff[0][i][j]) 
        
        #h_abij * a+_a * a+_b * a_i * a_j
        if self.coeff[1] is not None: 
            for i, j in [[_i, _j] for _j in range(self.n_qubits) for _i in range(_j)]:
                 for a, b in [[_a, _b] for _b in range(self.n_qubits) for _a in range(_b)]:
                     h += mp.two_body(self.n_qubits, a, b, i, j, coeff=self.coeff[1][a][b][i][j]) 
        h = h.reduce()
        pauli_list = []
        for pauli_sum_op in h:
            
            pauli_op = pauli_sum_op.to_pauli_op()
            if True:#str(pauli_op.primitive)=='XXI':
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
            
        
        #for p in pauli_list:
            #print(p.pauli_string)
        return pauli_list
            

class Pauli: #PauliOp with metods
    
    def __init__(self, pauli):
        self.pauli = pauli #PauliOp object
        self.pauli_string = str(pauli.primitive) #corresponding string
        self.coeff = pauli.coeff
        
    def pauli_to_qc(self, qc, qb):
        #basis change 
        pauli_string = self.pauli_string[::-1]
        for i, item in enumerate(pauli_string):
            if pauli_string[i] == 'X':
                qc.h(qb[i])
                
            elif pauli_string[i] == 'Y':
                qc.h(qb[i])
                qc.sdg(qb[i])
                
        return qc
    
    def apply_ham(self, qc, qb):
        #basis change 
        pauli_string = self.pauli_string[::-1]
        for i, item in enumerate(pauli_string):
            if pauli_string[i] == 'X':
                qc.x(qb[i])
                
            elif pauli_string[i] == 'Y':
                qc.y(qb[i])
                
            elif pauli_string[i] == 'Z':
                qc.z(qb[i])
                
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
            #print(state,value,sign * value/shots)
        return exp_value

    
    
    
    
    
    
    
    
    
    
    
    