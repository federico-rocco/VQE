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
from ansatz import Ansatz
from solver import Eigensolver
from algorithm import Algorithm
from optimizer import Optimizer
import numpy as np 
import matplotlib.pyplot as plt

orbitals = 3
fermions = 0
x=[]
y=[]
J = 1
xs = (x * 0.1 for x in range(0, 30))
for B in xs:


    ising_coeff = {'B':B,
                   'J':J
                   }
    hamiltonian = Hamiltonian(n_qubits=orbitals, coeff=ising_coeff)
    ansatz = Ansatz('ising', n_qubits=orbitals)
    
    options = {
        'shots':1024,
        'ibmq':False,
        'seed':10,
        'backend':'qasm_simulator', #qasm/aer/ibmq_something
        #'device':'ibm_nairobi' #to be simulated if using simulator
        }
    algorithm = Algorithm(options)
    optimizer = Optimizer('cobyla', disp=False, max_iter=1000)
    vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(ising=True), optimizer, algorithm)
    
    import time
    start_time = time.time()
    
    
    
    #ground state
    vqe.set_ising_params(ising_coeff)
    optimized_parameters = vqe.optimize_parameters(vqe.ising_expval)
    eigenvalue = vqe.ising_expval(optimized_parameters)/orbitals
    eigenstate = vqe.get_eigenstate(optimized_parameters)
    x.append(B)
    y.append(eigenvalue)

fig, ax = plt.subplots()

plt.ylabel('Energy')
plt.xlabel('B')
plt.scatter(x, y, marker=".", linewidth=2, zorder=3)
#print("result: ", eigenvalue, " expected: ", (-B-J), eigenstate)

#print("--- %s seconds ---" % (time.time() - start_time))

