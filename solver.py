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
import mapping as mp



        
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
        self.optimizer.set_loss_function(self.vqe_expval)      
        
        


    def vqe_qiskit(self):

        self.ansatz.count_excitations()
        qb = QuantumRegister(self.n_qubits)
        qc = QuantumCircuit(qb)
        theta = ParameterVector('Θ', self.ansatz.singles + self.ansatz.doubles)
        #self.ansatz.make_parameters()
        ansatz = self.ansatz(qc, qb, theta=theta)
        op = ~StateFn(self.hamiltonian)@StateFn(ansatz)
        grad = Gradient(grad_method='lin_comb').convert(operator=op)
        qi_sv = QuantumInstance(backend=self.algorithm.backend,
                        shots=self.algorithm.shots,
                        seed_simulator=2,
                        seed_transpiler=2)
        optimizer = SPSA(maxiter=100)  
        vqe = VQE(ansatz, optimizer=optimizer, gradient=grad, quantum_instance=qi_sv)
        return vqe.compute_minimum_eigenvalue(self.hamiltonian)
    

    def vqe_expval(self, theta=None):
        
        qb = QuantumRegister(self.n_qubits)
        cb = ClassicalRegister(self.n_qubits)
        qc = QuantumCircuit(qb, cb)
        #print(qc)
        
        ansatz_qc = self.ansatz(qc, qb, theta)
        
        
        E = 0
        for pauli_string in self.hamiltonian:
            qc = QuantumCircuit(qb, cb)
            pauli_op = pauli_string.to_pauli_op()
            qc = mp.pauli_to_qc(pauli_op, qc, qb)
            qc = qc + ansatz_qc
            #print(qc)
            #theta_val = 0*2*np.pi

            #circ = [qc.bind_parameters({self.ansatz.parameters: theta_val})
                   # for theta_val in theta_range]
            measurement = self.algorithm.measure(qc, qb, cb)
            expectation = self.algorithm.expectation(measurement)
            #print(pauli_op.coeff,expectation,pauli_op)
            E += pauli_op.coeff*expectation
        return E
    
    def optimize_parameters(self, theta=None):
        if theta is None:
            theta = self.ansatz.parameters
        params = self.optimizer(theta)
        return params
    
    
            
            



