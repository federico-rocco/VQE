# -*- coding: utf-8 -*-
"""
Created on Sat Jun 18 09:55:13 2022

@author: cosmo
"""
import numpy as np
import mapping as mp
from qiskit.circuit import ParameterVector
from qiskit import QuantumCircuit, QuantumRegister



class ansatz:
    
    def __init__(self, method, n_fermions=None, n_qubits=None, depth=1):
        self.n_fermions = n_fermions
        self.n_qubits = n_qubits
        self.singles = 0
        self.doubles = 0
        self.depth = depth
        self.method = method
        if self.method != 'UCCSD' and (self.n_fermions != 1 or self.n_qubits != 3):
            raise ValueError('This ansatz has a fixed number of fermions and orbitals')
        self.parameters = []
        self.mp2 = False
        
        
    def __call__(self,qc,qb,theta=None,qa=None,i=None):
        #builds and returns the ansatz circuit
        
        if theta is None:
            theta = self.new_parameters()
        
        if self.method == 'quarkonium':
            qc = self.quarkonium(qc, qb, theta)
        elif self.method == 'UCCSD':
            qc = self.UCCSD(qc, qb, theta)
        elif self.method == 'ising':
            qc = self.ising(qc, qb, theta)
        return qc
    
        
    
    def quarkonium(self, qc, qb, theta):
        #optimized ansatz from the quarkonium paper
        
        qc.ry(theta[1], 1)
        qc.ry(2*theta[0], 2)
        qc.cx(2, 0)
        qc.cx(0, 1)
        qc.x(2)
        qc.ry(-1*theta[1], 1)
        qc.cx(0, 1)
        qc.cx(1, 0)
        return qc 
    
    def ising(self, qc, qb, theta):
        
        qc.rz(theta[0], 0)
        qc.ry(theta[1], 0)

        # apply series of CNOT gates
        for i in range(1, self.n_qubits):
            qc.cx(0, i)

        # add parametrized single-qubit rotations around y
        for i in range(self.n_qubits):
            qc.ry(theta[2], i)
            
        return qc        
    
        
    
    def UCCSD(self, qc, qb, theta):
        #generic UCC ansatz

        qc = self.HartreeFock(qc, qb) #prepare HartreeFock: put a particle in the first orbital
        k = 0
        
        #single excitations
        for i in range(self.n_fermions):       
            for a in range(self.n_fermions, self.n_qubits):    
                qc = self.pauli_sum_op_to_exp_op_circuit(qc, qb, mp.one_body(self.n_qubits, a, i) - mp.one_body(self.n_qubits, i, a), theta[k])
                k += 1

        #double excitations
        for j in range(self.n_fermions):
                for i in range(j+1,self.n_fermions):
                    for b in range(self.fermions,self.n_qubits):
                        for a in range(b+1,self.n_qubits):
                            qc = self.pauli_sum_op_to_exp_op_circuit(qc, qb, mp.two_body(self.n_qubits, a, b, i, j) - mp.two_body(self.n_qubits, j, i, b, a), theta[k])
                            k += 1
        
        return qc
    
    def HartreeFock(self, qc, qb):
        for i in range(self.n_fermions):
            qc.x(qb[i])
        return qc
                       
    
    def pauli_sum_op_to_exp_op_circuit(self, qc, qb, pauli_sum_op, param):
        #from a PauliSumOp take each string (e.g. 'XYZ') and apply e^(param*string) to the circuit
        pauli_op_list = pauli_sum_op.to_pauli_op()
        for pauli_op in pauli_op_list:
            string = pauli_op.primitive.to_label()
            qc = self.exp_op(qc, qb, string, param)
        return qc
    
    
    def exp_op(self, qc, qb, pauli_string, theta): #exp(pauli_string)->circuit
    #e.g.: e^(theta*XYZ)
    #https://www.tcs.tifr.res.in/~pgdsen/pages/courses/2007/quantalgo/lectures/lec07.pdf
    #pag 30-31 https://www.researchgate.net/publication/233947759_The_Bravyi-Kitaev_transformation_for_quantum_computation_of_electronic_structure
        if not isinstance(pauli_string, str):
            raise ValueError("not a string")

        #rotate in the correct basis if needed
        for (i, letter) in zip(range(self.n_qubits), pauli_string):
            if letter == "X":
                qc.h(qb[i])
            elif letter == "Y":
                qc.h(qb[i])
                qc.rx(np.pi/2, qb[i])
    
            
        #compute parity and store it in the last qubit
        for i in range(self.n_qubits-1):
            qc.cx(qb[i], self.n_qubits-1)
        
        #apply parametrized rotation on the last qubit
        qc.rz(theta, qb[self.n_qubits-1])
        
        #uncompute parity
        for i in range(self.n_qubits-2, -1, -1):
            qc.cx(qb[i], qb[self.n_qubits-1])
    
        #restore Z basis
        for (i, letter) in zip(range(self.n_qubits-1, -1, -1), reversed(pauli_string)):
            if letter == "X":
                qc.h(qb[i])
            elif letter == "Y":
                qc.rx(-np.pi/2, qb[i])
                qc.h(qb[i])

        return qc
    
    
    def count_excitations(self): 
        #only depends on the number of orbitals and fermions
        
        for i in range(self.n_fermions):
            for a in range(self.n_fermions, self.n_qubits):
                self.singles += 1
        for j in range(self.n_fermions):
                for i in range(j+1,self.n_fermions):
                    for b in range(self.fermions,self.n_qubits):
                        for a in range(b+1,self.n_qubits):
                            self.doubles += 1
                            
                            
    def new_parameters(self,h=[0],e=[]):
        #generates random initial parameters (or mp2 parameters, to be implemented) 
    
        if self.method == "ising":
            return 2*np.pi*np.random.rand(3)
        
        if len(h)!=0 and len(e)!=0:
            self.mp2 = True
            
        params = []
        
        #t_i^a = 0
        for i in range(self.n_fermions):
            for a in range(self.n_fermions,self.n_qubits):
                self.singles += 1
                if self.mp2:
                    params.append(0)

        #compute t_ij^ab
        for i in range(self.n_fermions):
            for j in range(i+1,self.n_fermions):
                for a in range(self.n_fermions,self.n_qubits):
                    for b in range(a+1,self.n_qubits):
                        self.doubles += 1
                        if self.mp2:
                            t = (h[i,j,b,a]-h[i,j,a,b])/(e[i,i]+e[j,j]-e[a,a]-e[b,b])
                            params.append(t)

        if self.mp2:
            print('Using mp2 parameters')
            parameters = np.asarray(params)
        else:
            parameters = 2*np.pi*np.random.rand((self.singles+self.doubles)*self.depth)
        self.parameters = parameters
        return parameters
    
    

        

        

        
   
        
        
        
        
        
        
        
        
        
        
        