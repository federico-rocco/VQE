# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 18:49:15 2022

@author: cosmo
"""

from hamiltonian import Hamiltonian
from ansatz import UCC
from solver import Eigensolver, State
from qiskit import Aer, IBMQ
from algorithm import Algorithm
from optimizer import Minimizer
from quarkonium import *


coeffs = coeff()
hamiltonian = Hamiltonian(fermions,orbitals,coeffs)
ansatz = UCC(fermions,orbitals,method='quarkonium')

options = {
    'shots':1024,
    'ibmq':False,
    'backend':'aer_simulator_statevector', #qasm/aer/ibmq_something
    #'device':'ibm_nairobi' #to be simulated if using simulator
    }
algorithm = Algorithm(options)
optimizer = Minimizer('spsa', disp=False)

vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(), optimizer, algorithm)


#ground state
optimized_parameters = vqe.optimize_parameters(vqe.vqe_expval)
eigenvalue = vqe.vqe_expval(optimized_parameters)
eigenstate = vqe.get_eigenstate(optimized_parameters)
vqe.save_state(State(eigenvalue, eigenstate, optimized_parameters))
print("Ground state")
print("Energy: ", eigenvalue)
print("Eigenstate: ", eigenstate)
print("------------------------")



#first excited
optimized_parameters = vqe.optimize_parameters(vqe.expval_excited_state)
eigenvalue = vqe.vqe_expval(optimized_parameters)
eigenstate = vqe.get_eigenstate(optimized_parameters)
vqe.save_state(State(eigenvalue, eigenstate, optimized_parameters))
print("First excited state")
print("Energy: ", eigenvalue)
print("Eigenstate: ", eigenstate)
print("------------------------")



#second excited
optimized_parameters = vqe.optimize_parameters(vqe.expval_excited_state)
eigenvalue = vqe.vqe_expval(optimized_parameters)
eigenstate = vqe.get_eigenstate(optimized_parameters)
vqe.save_state(State(eigenvalue, eigenstate, optimized_parameters))
print("Second excited state state")
print("Energy: ", eigenvalue)
print("Eigenstate: ", eigenstate)
print("------------------------")