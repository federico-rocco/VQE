# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 09:55:13 2022

@author: cosmo
"""
import numpy as np
import mapping as mp
from qiskit.circuit import ParameterVector
from qiskit import QuantumCircuit, QuantumRegister



class UCC:
    
    def __init__(self,f,q,depth=1,trunc='UCCSD',quark=False):
        self.n_fermions = f
        self.n_qubits = q
        self.singles = 0
        self.doubles = 0
        self.depth = depth
        self.trunc = trunc
        self.parameters = []
        self.quarkonium = quark
        self.mp2 = False
        
        
    def __call__(self,qc,qb,theta=None,qa=None,i=None):
        if theta is None:
            theta = self.new_parameters()
            print("Initial parameters: ", theta)
        qc = self.UCCSD(qc, qb, theta)
        return qc
                             
    
    def pauli_sum_op_to_exp_op_circuit(self, qc, qb, qa, pauli_sum_op, param):
        pauli_op_list = pauli_sum_op.to_pauli_op()
        for pauli_op in pauli_op_list:
            string = pauli_op.primitive.to_label()
            qc = self.exp_op(qc, qb, qa, string, param)
        return qc
    
    
    def exp_op(self, qc, qb, qa, pauli_string, theta): #exp(pauli_string)->circuit
    #e.g.: e^(theta*XYZ)
    #https://www.tcs.tifr.res.in/~pgdsen/pages/courses/2007/quantalgo/lectures/lec07.pdf
    #pag 30-31 https://www.researchgate.net/publication/233947759_The_Bravyi-Kitaev_transformation_for_quantum_computation_of_electronic_structure
        if not isinstance(pauli_string, str):
            raise ValueError("not a string")

        #rotate in the correct basis if needed
        for (i, letter) in zip(range(len(qb)-1), pauli_string):
            if letter == "X":
                qc.h(qb[i])
            elif letter == "Y":
                qc.rx(-np.pi/2, qb[i])
        
        
        #compute parity with CNOTs and store it in the ancilla
        for qubit in qb:
            qc.cx(qubit, qa)
        
        #apply parametrized Rz rotation to the last qubit
        qc.rz(theta, qa)
        
        #uncompute parity
        for qubit in reversed(qb):
            qc.cx(qubit, qa)
    
        #restore Z basis
        for (i, letter) in zip(range(len(qb)-1, -1, -1), reversed(pauli_string)):
            if letter == "X":
                qc.h(qb[i])
            elif letter == "Y":
                qc.rx(np.pi/2, qb[i])

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
                            
                            
    def new_parameters(self,h=[],e=[]):
    
        if len(h)!=0 and len(e)!=0:
            self.mp2 = True
            
        params = []
        if 'S' in self.trunc: #t_i^a = 0
            for i in range(self.n_fermions):
                for a in range(self.n_fermions,self.n_qubits):
                    self.singles += 1
                    if self.mp2:
                        params.append(0)

        if 'D' in self.trunc: #compute t_ij^ab
            for i in range(self.n_fermions):
                for j in range(i+1,self.n_fermions):
                    for a in range(self.n_fermions,self.n_qubits):
                        for b in range(a+1,self.n_qubits):
                            self.doubles += 1
                            if self.mp2:
                                t = (h[i,j,b,a]-h[i,j,a,b])/(e[i,i]+e[j,j]-e[a,a]-e[b,b])
                                params.append(t)

        if self.mp2:
            parameters = np.asarray(params)
        else:
            parameters = 2*np.pi*np.random.rand((self.singles+self.doubles)*self.depth)
        self.parameters = parameters
        return parameters#ParameterVector('Θ', self.singles + self.doubles)
    
    
    def UCCSD(self, qc, qb, theta):

        if self.quarkonium:

            qc.ry(theta[1], 1)
            qc.ry(2*theta[0], 2)
            qc.cx(2, 0)
            qc.cx(0, 1)
            qc.x(2)
            qc.ry(-1*theta[1], 1)
            qc.cx(0, 1)
            qc.cx(1, 0)
            return qc 
        
        else:
            
            
            #[qc.x(qb[i]) for i in range(self.n_qubits)]
            qc.x(qb[0])
            #qc.x(qb[1])
            #qc.x(qb[2])
            qa = QuantumRegister(1)
            qc.add_register(qa)
            k = 0
            
            for i in range(self.n_fermions):       
                for a in range(self.n_fermions, self.n_qubits):
                    
                    qc = self.pauli_sum_op_to_exp_op_circuit(qc, qb, qa, mp.one_body(self.n_qubits, a, i) - mp.one_body(self.n_qubits, i, a), theta[k])
                    k += 1
    
            for j in range(self.n_fermions):
                    for i in range(j+1,self.n_fermions):
                        for b in range(self.fermions,self.n_qubits):
                            for a in range(b+1,self.n_qubits):
                                qc = self.pauli_sum_op_to_exp_op_circuit(qc, qb, qa, mp.two_body(self.n_qubits, a, b, i, j) - mp.two_body(self.n_qubits, j, i, b, a), theta[k])
                                k += 1
        
            #print(qc)
            return qc
        

        

        
   
        
        
        
        
        
        
        
        
        
        
        