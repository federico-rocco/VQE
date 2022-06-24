# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 09:55:13 2022

@author: cosmo
"""
import numpy as np
import mapping as mp

from qiskit.circuit import QuantumCircuit, ParameterVector

class UCC:
    
    def __init__(self,n,l):
        self.n_fermions = n
        self.n_qubits = l
        self.singles = 0
        self.doubles = 0
        self.parameters = []
                             
    
    def pauli_sum_op_to_exp_op_circuit(self, pauli_sum_op, param):
        qc = QuantumCircuit(self.n_qubits)
        pauli_op_list = pauli_sum_op.to_pauli_op()
        for pauli_op in pauli_op_list:
            string = pauli_op.primitive.to_label()
            qc += self.exp_op(string, param)
        return qc
    
    
    def exp_op(self, pauli_string, theta): #exp(pauli_string)->circuit
        if not isinstance(pauli_string, str):
            raise ValueError("not a string")
        qc = QuantumCircuit(self.n_qubits)

        for (qubit, letter) in zip(range(self.n_qubits), pauli_string):
            if letter == "X":
                qc.h(qubit)
            elif letter == "X":
                qc.rx(-np.pi/2, qubit)
                
        for qubit in range(self.n_qubits-1):
            qc.cx(qubit, self.n_qubits-1)
            
        qc.rz(theta, self.n_qubits-1)
        
        for qubit in range(self.n_qubits-2, -1, -1):
            qc.cx(qubit, self.n_qubits-1)
    
        for (qubit, letter) in zip(range(self.n_qubits-1, -1, -1), reversed(pauli_string)):
            if letter == "X":
                qc.h(qubit)
            elif letter == "Y":
                qc.rx(np.pi/2, qubit)

        return qc
    
    
    def count_excitations(self):
        
        for i in range(self.n_fermions):
            for a in range(self.n_fermions, self.n_qubits):
                self.singles += 1
        for j in range(self.n_fermions):
                for i in range(j+1,self.n_fermions):
                    for b in range(self.fermions,self.n_qubits):
                        for a in range(b+1,self.n_qubits):
                            self.doubles += 1
                            
                            
    def make_parameters(self):                        
        self.parameters = ParameterVector('theta', self.singles+self.doubles)
    
    def UCC_ansatz(self):
        self.count_excitations()
        self.make_parameters()
        
        ansatz = QuantumCircuit(self.n_qubits)
        k = 0
        
        for i in range(self.n_fermions):       
            for a in range(self.n_fermions, self.n_qubits):
                ansatz += self.pauli_sum_op_to_exp_op_circuit(mp.one_body(self.n_qubits, a, i) - mp.one_body(self.n_qubits, i, a), self.parameters[k])
                k += 1

        for j in range(self.n_fermions):
                for i in range(j+1,self.n_fermions):
                    for b in range(self.fermions,self.n_qubits):
                        for a in range(b+1,self.n_qubits):
                                ansatz += self.pauli_sum_op_to_exp_op_circuit(mp.two_body(self.n_qubits, a, b, i, j) - mp.two_body(self.n_qubits, j, i, b, a), self.parameters[k])
    
        return ansatz
 
      
        
        
        
        
        
        
        
        
        
        
        