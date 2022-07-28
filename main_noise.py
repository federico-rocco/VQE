# -*- coding: utf-8 -*-
"""
Created on Thu May 12 10:59:48 2022

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
    'ibmq':True,
    'backend':'ibm_nairobi', #qasm_simulator/aer_simulator_statevector/ibmq_something
    #'device':'ibm_nairobi' #to be simulated if using simulator
    }
algorithm = Algorithm(options)
optimizer = Minimizer('spsa', max_iter=50)

vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(), optimizer, algorithm)


#ground state
print("no folding, normal noise")
#vqe.set_folding(0)
optimized_parameters = vqe.optimize_parameters(vqe.vqe_expval)
print("Optimized: ", vqe.vqe_expval(optimized_parameters))

print("------------------------")



print("n_folding = 1")
vqe.set_folding(1)
optimized_parameters = vqe.optimize_parameters(vqe.vqe_expval)
print("Optimized: ", vqe.vqe_expval(optimized_parameters))

print("------------------------")



print("n_folding = 2")
vqe.set_folding(2)
optimized_parameters = vqe.optimize_parameters(vqe.vqe_expval)
print("Optimized: ", vqe.vqe_expval(optimized_parameters))

print("------------------------")



print("n_folding = 3")
vqe.set_folding(3)
optimized_parameters = vqe.optimize_parameters(vqe.vqe_expval)
print("Optimized: ", vqe.vqe_expval(optimized_parameters))

print("------------------------")
