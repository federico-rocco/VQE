# -*- coding: utf-8 -*-
"""
Created on Mon Aug 22 15:59:10 2022

@author: cosmo
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 19:22:37 2022

@author: cosmo
"""

from hamiltonian import Hamiltonian
from ansatz import ansatz
from solver import Eigensolver
from algorithm import Algorithm
from optimizer import Optimizer
import numpy as np 


orbitals = 3
fermions = 0
B = 0.1
J = 1
ising_coeff = {'B':B,
               'J':J
               }
hamiltonian = Hamiltonian(n_qubits=orbitals, coeff=ising_coeff)
ansatz = ansatz('ising', n_qubits=orbitals)

options = {
    'shots':1024,
    'ibmq':False,
    'seed':10,
    'backend':'aer_simulator_statevector', #qasm/aer/ibmq_something
    #'device':'ibm_nairobi' #to be simulated if using simulator
    }
algorithm = Algorithm(options)
optimizer = Optimizer('spsa', disp=False, max_iter=100)
vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(ising=True), optimizer, algorithm)

import time
start_time = time.time()



#ground state
optimized_parameters = vqe.optimize_parameters(vqe.vqe_expval)
eigenvalue = vqe.vqe_expval(optimized_parameters)
eigenstate = vqe.get_eigenstate(optimized_parameters)
print("result: ", eigenvalue, " expected: ", (-B/2-J/4)*orbitals, eigenstate)

print("--- %s seconds ---" % (time.time() - start_time))