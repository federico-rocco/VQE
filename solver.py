# -*- coding: utf-8 -*-
"""
Created on Thu May 12 11:45:00 2022

@author: cosmo
"""


from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit.algorithms import VQE
from qiskit.algorithms.optimizers import SPSA #COBYLA, CG, SLSQP 
from qiskit.utils import QuantumInstance
from qiskit.opflow import StateFn
from qiskit.opflow.gradients import Gradient
from qiskit.circuit import ParameterVector, Parameter
import numpy as np
import copy



        
class Eigensolver:
    
    def __init__(self, 
                 n_fermions, 
                 n_qubits, 
                 ansatz, 
                 hamiltonian, 
                 optimizer, 
                 algorithm):    
        self.n_fermions = n_fermions
        self.n_qubits = n_qubits
        self.ansatz = ansatz
        self.hamiltonian = hamiltonian
        self.algorithm = algorithm
        self.optimizer = optimizer
        self.i = 0
        self.solved_states = []
        self.n_folding = 0
        
        


    def vqe_qiskit(self):     
        
        
        self.ansatz.count_excitations()
        qb = QuantumRegister(self.n_qubits)
        qc = QuantumCircuit(qb)
        theta = ParameterVector('Θ', self.ansatz.singles + self.ansatz.doubles)
        ansatz = self.ansatz(qc, qb, theta=theta)
        
        op = 0
        ham = 0
        for pauli in self.hamiltonian:
            op += ~StateFn(pauli.pauli.to_pauli_op())@StateFn(ansatz)
            ham += pauli.pauli
            
        grad = Gradient(grad_method='lin_comb').convert(operator=op, params=theta)
        qi_sv = QuantumInstance(backend=self.algorithm.backend,
                        shots=self.algorithm.shots,
                        seed_simulator=2,
                        seed_transpiler=2) 
        vqe = VQE(ansatz, optimizer=self.optimizer.opt, gradient=grad, quantum_instance=qi_sv)
        return vqe.compute_minimum_eigenvalue(ham)
    

    def vqe_expval(self, theta=None):
        qb = QuantumRegister(self.n_qubits)
        cb = ClassicalRegister(self.n_qubits)
        qc = QuantumCircuit(qb)
        ansatz_qc = self.ansatz(qc, qb, theta)

        E = 0
        for pauli in self.hamiltonian:
            qc = QuantumCircuit(qb, cb)
            qc = ansatz_qc + qc
            qc = pauli.pauli_to_qc(qc, qb)
            qc = self.fold(qc)
            measurement = self.algorithm.measure(qc, qb, cb)
            expectation = pauli.expectation(measurement, self.algorithm.shots)
            E += pauli.coeff*expectation

        return E
    
    def fold(self, qc):
        
        qc0 = copy.deepcopy(qc)
        qc_inv = qc.inverse()
        for i in range(self.n_folding):
            qc += (qc_inv + qc0)   
        return qc
    
    def set_folding(self, n):
        self.n_folding = n
    
    def optimize_parameters(self, loss_function, theta=None):
        if theta is None:
            theta = self.ansatz.new_parameters()

        result = self.optimizer(loss_function, theta)
        return result
    
    def measure_ansatz(self, theta=None):
        qb = QuantumRegister(self.n_qubits)
        cb = ClassicalRegister(self.n_qubits)
        qc = QuantumCircuit(qb, cb)
        qc = self.ansatz(qc, qb, theta)
        measurement = self.algorithm.measure(qc, qb, cb)
        return measurement
            
    def get_eigenstate(self, theta=None) :       
        qb = QuantumRegister(self.n_qubits)
        cb = ClassicalRegister(self.n_qubits)
        qc = QuantumCircuit(qb, cb)

        ansatz_qc = self.ansatz(qc, qb, theta)
        measurement = self.algorithm.measure(ansatz_qc, qb, cb)
        return measurement


    def save_state(self, state):
        self.solved_states.append(state)
        
        
    def find_delta(self):
        delta = 0
        for pauli in self.hamiltonian:
            delta += 2*abs(pauli.coeff)
        return delta
   
    def expval_excited_state(self, theta_k):
        
        qb = QuantumRegister(self.n_qubits)
        cb = ClassicalRegister(self.n_qubits)
        qc = QuantumCircuit(qb, cb)
        ansatz_k = self.ansatz(qc, qb, theta_k)
        overlap = 0
        for state in self.solved_states:
            qc = QuantumCircuit(qb, cb)
            inverse = self.ansatz(qc, qb, state.parameters).inverse()
            
            new_qc = ansatz_k + inverse
            
            measurement = self.algorithm.measure(new_qc, qb, cb)

            zero = '0'
            for q in range(self.n_qubits-1):
                zero += '0'
            if zero in measurement:
                overlap += measurement[zero]

        return self.find_delta()*overlap + self.vqe_expval(theta_k)


        
class State:

    def __init__(self, eigenvalue, eigenstate, parameters):
        self.eigenvalue = eigenvalue
        self.eigenstate = eigenstate
        self.parameters = parameters
        
    