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
from folding import circuit_folding, gate_folding
import numpy as np
import copy
from qiskit.opflow import PauliExpectation


        
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
        self.solved_states = []
        self.n_folding = 1
        self.factor = 1
        
        


    def vqe_qiskit(self):     
        #uses qiskit.algorithms.VQE 
        
        self.ansatz.count_excitations()
        qb = QuantumRegister(self.n_qubits)
        qc = QuantumCircuit(qb)
        theta = ParameterVector('Θ', self.ansatz.singles + self.ansatz.doubles) #qiskit requires a parametrized ansatz
        ansatz = self.ansatz(qc, qb, theta=theta)
        ansatz = circuit_folding(ansatz, scaling=self.n_folding)
        
        #operator = sum <H|psi>
        operator = 0
        hamiltonian = 0
        for pauli in self.hamiltonian:
            operator += ~StateFn(pauli.pauli.to_pauli_op())@StateFn(ansatz)
            hamiltonian += pauli.pauli
            
        grad = Gradient(grad_method='lin_comb').convert(operator=operator, params=theta)
        qi_sv = QuantumInstance(backend=self.algorithm.backend,
                        shots=self.algorithm.shots) 
        vqe = VQE(ansatz, optimizer=self.optimizer.opt, quantum_instance=qi_sv, gradient=grad)
        return vqe.compute_minimum_eigenvalue(hamiltonian)
    

    def expval(self, theta=None):
        #self-made eigensolver, not using qiskit
        #build ansatz, then add a term of the hamiltonian at a time and compute the
        #expectation value <psi|H|psi>
        
        qb = QuantumRegister(self.n_qubits)
        cb = ClassicalRegister(self.n_qubits)
        qc = QuantumCircuit(qb)
        ansatz_qc = self.ansatz(qc, qb, theta)
        ansatz_qc = circuit_folding(ansatz_qc, scaling=self.n_folding) #extend the circuit in case of noise extrapolation

        E = 0
        for pauli in self.hamiltonian:
            qc = QuantumCircuit(qb, cb)
            qc = qc.compose(ansatz_qc) 
            qc = pauli.pauli_to_qc(qc, qb) #the piece of hamiltonian acts on psi
            measurement = self.algorithm.measure(qc, qb, cb) #measure the circuit
            expectation = pauli.expectation(measurement, self.algorithm.shots) #using relative frequencies compute exp_value
            E += pauli.coeff*expectation

        return E
    
    def analytic_expval(self, theta=None):
        #build ansatz, then add a term of the hamiltonian at a time and compute the
        #expectation value <psi|H|psi> analytically using qiskit tools
        
        qb = QuantumRegister(self.n_qubits)
        qc = QuantumCircuit(qb)
        ansatz_qc = self.ansatz(qc, qb, theta)
        ansatz_qc = circuit_folding(ansatz_qc, scaling=self.n_folding) #extend the circuit in case of noise extrapolation

        E = 0
        for pauli in self.hamiltonian:
            expectation = ~StateFn(pauli.pauli.to_pauli_op())@StateFn(ansatz_qc)
            E += expectation.eval().real

        return E
        
    
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
        self.n_folding=5
        qc=self.fold(qc)
        measurement = self.algorithm.measure(qc, qb, cb)
        return measurement
            
    def get_eigenstate(self, theta=None) :       
        qb = QuantumRegister(self.n_qubits)
        cb = ClassicalRegister(self.n_qubits)
        qc = QuantumCircuit(qb, cb)

        ansatz_qc = self.ansatz(qc, qb, theta)
        measurement = self.algorithm.measure(ansatz_qc, qb, cb)
        return measurement
    
    def get_ansatz(self, theta=None) :       
        qb = QuantumRegister(self.n_qubits)
        qc = QuantumCircuit(qb)

        ansatz_qc = self.ansatz(qc, qb, theta)
        return ansatz_qc


    def save_state(self, state):
        self.solved_states.append(state)
        
        
    def find_delta(self):
        delta = 0 #8715
        for pauli in self.hamiltonian:
            delta += 2*abs(pauli.coeff)
        return delta/self.factor
   
    def expval_excited_state(self, theta_k):
        #https://arxiv.org/pdf/1805.08138.pdf
        #F(Θk) = <ψ(Θk)| H |ψ(Θk)> + sum_i βi * |<ψ(Θk)|ψ(Θi)>|^2
        #the first term is the normal expectation value, like for the ground state
        #the second term is the overlap with all the previous excited states 
        #it has to be minimized since the state k should be orthogonal to all the others
        #βi can be set = to a very big value delta for all the states

        qb = QuantumRegister(self.n_qubits)
        cb = ClassicalRegister(self.n_qubits)
        qc = QuantumCircuit(qb, cb)
        ansatz_k = self.ansatz(qc, qb, theta_k)
        
        overlap = 0
        for state in self.solved_states:
            qc = QuantumCircuit(qb, cb)
            inverse = self.ansatz(qc, qb, state.parameters).inverse()
            
            new_qc = ansatz_k + inverse #<ψ(Θi)|ψ(Θk)> 
            
            measurement = self.algorithm.measure(new_qc, qb, cb)

            zero = '0'*self.n_qubits
            if '0'*self.n_qubits in measurement:
                overlap += measurement[zero] 
                #|<ψ(Θi)|ψ(Θk)>|^2 = |<0|ψ(Θi)^-1 ψ(Θk)|0>|^2
                #so we just count the all '0' qubits

        return self.find_delta()*overlap + self.expval(theta_k)


        
class State:

    def __init__(self, eigenvalue, eigenstate, parameters):
        self.eigenvalue = eigenvalue
        self.eigenstate = eigenstate
        self.parameters = parameters
        
    