# -*- coding: utf-8 -*-
"""
Created on Thu May 12 10:59:48 2022

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


coeffs = coeff()
hamiltonian = Hamiltonian(fermions,orbitals,coeffs)
ansatz = ansatz('quarkonium', n_fermions=fermions, n_qubits=orbitals)


options = {
    'seed':1,
    'shots':1024,
    'ibmq':False,
    'backend':'qasm_simulator', #qasm/aer/ibmq_something
    'device':'ibm_nairobi' #to be simulated if using simulator
    }
algorithm = Algorithm(options)
optimizer = Optimizer('spsa', max_iter=1000)
vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(), optimizer, algorithm)



import time
start_time = time.time()
executions=10

x = []
y = []
y_std= []
weights = []
for lamda in range(1,5,1):
    x.append(lamda)
    print("lamda = ", lamda)
    vqe.set_folding(lamda)
    energies = []
    for i in range(executions):
        optimized_parameters = vqe.optimize_parameters(vqe.vqe_expval)
        energy = vqe.vqe_expval(optimized_parameters)
        if energy < 700:
            energies.append(energy)
    mean_energy = np.mean(energies)
    std = np.std(energies)
    y.append(mean_energy)
    y_std.append(std)
    weights.append(1/std)
    print("Optimized energy: ", mean_energy, "all energies: ", energies)
    
    print("------------------------")

b, m = np.polynomial.polynomial.polyfit(x, y, 1, rcond=None, full=False)
plt.errorbar(x,y,yerr=y_std)
plt.xlim([-1,11])
plt.ylim([450,700])
x = np.array(x)
plt.plot(x, m*x + b)
plt.show()


print("--- %s seconds ---" % (time.time() - start_time))