# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 18:49:15 2022

@author: cosmo
"""

from hamiltonian import Hamiltonian
from ansatz import ansatz
from solver import Eigensolver, State
from qiskit import Aer, IBMQ
from algorithm import Algorithm
from optimizer import Optimizer
from quarkonium import *
import matplotlib.pyplot as plt


coeffs = coeff()
hamiltonian = Hamiltonian(fermions,orbitals,coeffs)
ansatz = ansatz('quarkonium', n_fermions=fermions, n_qubits=orbitals)

options = {
    'shots':1024,
    'ibmq':False,
    'backend':'qasm_simulator', #qasm/aer/ibmq_something
    #'device':'ibm_nairobi' #to be simulated if using simulator
    }
algorithm = Algorithm(options)
optimizer = Optimizer('spsa', disp=False)



x = []
y = []
for delta in range(600,1500,20):
    x.append(delta)
    wrong = 0
    wrong_eig = []
    for i in range(20):       
        vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(), optimizer, algorithm)
        vqe.factor = delta 
        #ground state
        optimized_parameters = vqe.optimize_parameters(vqe.vqe_expval)
        eigenvalue = vqe.vqe_expval(optimized_parameters)
        eigenstate = vqe.get_eigenstate(optimized_parameters)
        vqe.save_state(State(eigenvalue, eigenstate, optimized_parameters))
        
        
        
        #first excited
        optimized_parameters1 = vqe.optimize_parameters(vqe.expval_excited_state)
        eigenvalue1 = vqe.vqe_expval(optimized_parameters1)
        eigenstate1 = vqe.get_eigenstate(optimized_parameters1)
        vqe.save_state(State(eigenvalue1, eigenstate1, optimized_parameters1))

        
        
        
        #second excited
        optimized_parameters2 = vqe.optimize_parameters(vqe.expval_excited_state)
        eigenvalue2 = vqe.vqe_expval(optimized_parameters2)
        eigenstate2 = vqe.get_eigenstate(optimized_parameters2)
        vqe.save_state(State(eigenvalue2, eigenstate2, optimized_parameters2))

        
        if not eigenvalue2 > eigenvalue1 > eigenvalue:
            wrong += 1
            wrong_eig.append(eigenvalue)
            wrong_eig.append(eigenvalue1)
            wrong_eig.append(eigenvalue2)
            wrong_eig.append("---")

     
    y.append(wrong*5)       
    print("With delta/", delta, wrong*5, "% error.")
    print("========================")
plt.plot(x,y)
plt.show