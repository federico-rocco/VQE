# -*- coding: utf-8 -*-
"""
Created on Sat Oct 15 17:33:26 2022

@author: cosmo
"""

from hamiltonian import Hamiltonian
from ansatz import Ansatz
from solver import Eigensolver, State
from qiskit import Aer, IBMQ
from algorithm import Algorithm
from optimizer import Optimizer
from quarkonium import *
import numpy as np
import matplotlib.pyplot as plt
import os

import time
start_time = time.time()


fermions = 1
orbitals = 3
model = Quarkonium('charmonium')
coeffs = model.coeff(orbitals)
hamiltonian = Hamiltonian(fermions,orbitals,coeffs)
ansatz = Ansatz('quarkonium', n_fermions=fermions, n_qubits=orbitals)


options = {
    #'seed':1,
    'shots':1024*100,
    'ibmq':False,
    'backend':'qasm_simulator', #qasm/aer/ibmq_something
    'device':'ibmq_athens' #to be simulated if using simulator
    }
algorithm = Algorithm(options)
optimizer = Optimizer('spsa', max_iter=200)
vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(), optimizer, algorithm)
optimized_parameters = vqe.optimize_parameters(vqe.expval)

print('found parameters')

executions=100

x = []
y = []
y_std= []
weights = []
for lamda in range(1,6,1):
    x.append(lamda)
    print("lamda = ", lamda)
    vqe.set_folding(lamda)
    energies = []
    for i in range(executions):
        
        energy = vqe.expval(optimized_parameters)
        energies.append(energy)
    mean_energy = np.mean(energies)
    std = np.std(energies)
    y.append(mean_energy)
    y_std.append(std)
    #weights.append(1/std)
    print("Optimized energy: ", mean_energy)
    
    print("------------------------")

x = np.array(x)
y = np.array(y)
b, m = np.polynomial.polynomial.polyfit(x, y, 1, rcond=None, full=False)
#plt.errorbar(x,y,yerr=y_std)
#plt.xlim([-1,11])
#plt.ylim([450,700])

x2 = np.insert(x, 0, 0)

y2 = np.insert(y, 0, b)
y_std = np.array(y_std)
y_std2 = np.insert(y_std, 0, np.mean(y_std))

plt.errorbar(x,y,yerr=y_std, ecolor='tab:red', capsize=3, fmt="r--o",linewidth=1, mfc='none')
plt.plot(x2, y2, linestyle='dashed', color='red')
plt.plot([0],[b],'r*')
plt.fill_between(x2, y2-y_std2, y2+y_std2, edgecolor='pink', facecolor='moccasin')
plt.ylabel('E [MeV]')
plt.xlabel('λ')
plt.show()


print("--- %s seconds ---" % (time.time() - start_time))