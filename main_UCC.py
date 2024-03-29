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


coeffs = coeff()
hamiltonian = Hamiltonian(fermions,orbitals,coeffs)
ansatz = ansatz('UCCSD', n_fermions=fermions, n_qubits=orbitals)

options = {
    'shots':1024,
    'ibmq':False,
    'seed':10,
    'backend':'aer_simulator_statevector', #qasm/aer/ibmq_something
    #'device':'ibm_nairobi' #to be simulated if using simulator
    }
algorithm = Algorithm(options)
optimizer = Optimizer('spsa', disp=False, max_iter=10000)

vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(), optimizer, algorithm)

import time
start_time = time.time()



#ground state
optimized_parameters = vqe.optimize_parameters(vqe.vqe_expval)
eigenvalue = vqe.vqe_expval(optimized_parameters)
print(eigenvalue)

print("--- %s seconds ---" % (time.time() - start_time))