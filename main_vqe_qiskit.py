# -*- coding: utf-8 -*-
"""
Created on Fri Sep  2 17:30:47 2022

@author: cosmo
"""

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
import copy 


coeffs = coeff()
hamiltonian = Hamiltonian(fermions,orbitals,coeffs)
ansatz = ansatz('quarkonium', n_fermions=fermions, n_qubits=orbitals)

options = {
    'shots':1024,
    'ibmq':False,
    'seed':10,
    'backend':'qasm_simulator', #qasm/aer/ibmq_something
    'device':'ibmq_athens' #to be simulated if using simulator
    }
algorithm = Algorithm(options)
optimizer = Optimizer('spsa', disp=False, max_iter=50)

import time
start_time = time.time()


for n_folding in range(5):
    
    vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(), optimizer, algorithm)
    vqe.set_folding(n_folding)
    #ground state
    #optimized_parameters = vqe.optimize_parameters(vqe.vqe_expval)
    eigenvalue = vqe.vqe_qiskit().optimal_value
    print("folding: ", n_folding, "result:", eigenvalue )
    

#mpt.pyplot.hist(eigs)
print("--- %s seconds ---" % (time.time() - start_time))