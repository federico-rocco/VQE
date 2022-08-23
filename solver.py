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
        self.solved_states = []
        self.n_folding = 0
        
        


    def vqe_qiskit(self):     
        #uses qiskit.algorithms.VQE 
        
        self.ansatz.count_excitations()
        qb = QuantumRegister(self.n_qubits)
        qc = QuantumCircuit(qb)
        theta = ParameterVector('Θ', self.ansatz.singles + self.ansatz.doubles) #qiskit requires a parametrized ansatz
        ansatz = self.ansatz(qc, qb, theta=theta)
        
        #operator = sum <H|psi>
        operator = 0
        hamiltonian = 0
        for pauli in self.hamiltonian:
            operator += ~StateFn(pauli.pauli.to_pauli_op())@StateFn(ansatz)
            hamiltonian += pauli.pauli
            
        grad = Gradient(grad_method='lin_comb').convert(operator=operator, params=theta)
        qi_sv = QuantumInstance(backend=self.algorithm.backend,
                        shots=self.algorithm.shots,
                        seed_simulator=2,
                        seed_transpiler=2) 
        vqe = VQE(ansatz, optimizer=self.optimizer.opt, gradient=grad, quantum_instance=qi_sv)
        return vqe.compute_minimum_eigenvalue(hamiltonian)
    

    def vqe_expval(self, theta=None):
        #self-made eigensolver, not using qiskit
        #build ansatz, then add a term of the hamiltonian at a time and compute the
        #expectation value <psi|H|psi>
        
        qb = QuantumRegister(self.n_qubits)
        cb = ClassicalRegister(self.n_qubits)
        qc = QuantumCircuit(qb)
        ansatz_qc = self.ansatz(qc, qb, theta)

        E = 0
        for pauli in self.hamiltonian:
            qc = QuantumCircuit(qb, cb)
            qc = ansatz_qc + qc
            qc = pauli.pauli_to_qc(qc, qb) #the piece of hamiltonian acts on psi
            qc = self.fold(qc) #extend the circuit in case of noise extrapolation
            measurement = self.algorithm.measure(qc, qb, cb) #measure the circuit
            expectation = pauli.expectation(measurement, self.algorithm.shots) #using relative frequencies compute exp_value
            E += pauli.coeff*expectation

        return E
    
    def fold(self, qc):
        #U ---> U*(Udagger * U)^n where U is a unitary circuit
        
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
        return delta/1000
   
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

            zero = '0'
            for q in range(self.n_qubits-1):
                zero += '0'
            if zero in measurement:
                overlap += measurement[zero] 
                #|<ψ(Θi)|ψ(Θk)>|^2 = |<0|ψ(Θi)^-1 ψ(Θk)|0>|^2
                #so we just count the all '0' qubits

        return self.find_delta()*overlap + self.vqe_expval(theta_k)


        
class State:

    def __init__(self, eigenvalue, eigenstate, parameters):
        self.eigenvalue = eigenvalue
        self.eigenstate = eigenstate
        self.parameters = parameters
        
    