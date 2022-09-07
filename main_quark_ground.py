# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 19:22:37 2022

@author: cosmo
"""

from hamiltonian import Hamiltonian
from ansatz import ansatz
from solver import Eigensolver, State
from qiskit import Aer, IBMQ
from algorithm import Algorithm
from optimizer import Optimizer
from quarkonium import *
import matplotlib as mpt

coeffs = coeff()
hamiltonian = Hamiltonian(fermions,orbitals,coeffs)
ansatz = ansatz('quarkonium', n_fermions=fermions, n_qubits=orbitals)

options = {
    'shots':1024,
    'ibmq':True,
    'seed':10,
    'backend':'ibm_oslo', #qasm/aer/ibmq_something
    #'device':'ibm_nairobi' #to be simulated if using simulator
    }
algorithm = Algorithm(options)
optimizer = Optimizer('', disp=False, max_iter=10)

vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(), optimizer, algorithm)

import time
start_time = time.time()

eigs = []
for i in range(1):
    #ground state
    optimized_parameters = vqe.optimize_parameters(vqe.vqe_expval)
    eigenvalue = vqe.vqe_expval(optimized_parameters)
    print("singlet:", eigenvalue + VSS(0), "triplet:", eigenvalue + VSS(1))
    eigs.append(eigenvalue)

#mpt.pyplot.hist(eigs)
print("--- %s seconds ---" % (time.time() - start_time))