# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 15:44:42 2022

@author: cosmo
"""

from hamiltonian import Hamiltonian
from ansatz import ansatz
from solver import Eigensolver, State
from qiskit import Aer, IBMQ
from algorithm import Algorithm
from optimizer import Optimizer
from quarkonium import *
import numpy as np
import matplotlib.pyplot as plt
import os

coeffs = coeff()
hamiltonian = Hamiltonian(fermions,orbitals,coeffs)
ansatz = ansatz('quarkonium', n_fermions=fermions, n_qubits=orbitals)


options = {
    'seed':1,
    'shots':1024,
    'ibmq':False,
    'backend':'qasm_simulator', #qasm/aer/ibmq_something
    'device':'ibmq_athens' #to be simulated if using simulator
    }
algorithm = Algorithm(options)
optimizer = Optimizer('spsa', max_iter=100)
vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(), optimizer, algorithm)
optimized_parameters = vqe.optimize_parameters(vqe.vqe_expval)


import time
start_time = time.time()
x = []
y = []

for lamda in range(1,6,1):
    x.append(lamda)
    print("lamda = ", lamda)
    vqe.set_folding(lamda)
    energy = vqe.vqe_expval(optimized_parameters)
    y.append(energy)
    print("Optimized energy: ", energy)
    
    print("------------------------")

b, m = np.polynomial.polynomial.polyfit(x, y, 1, rcond=None, full=False)
x = np.array(x)
plt.plot(x, m*x + b)
plt.scatter(x, y)
plt.show()


print("--- %s seconds ---" % (time.time() - start_time))