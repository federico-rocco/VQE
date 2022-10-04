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
import matplotlib.pyplot as plt

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
    'backend':'qasm_simulator', #qasm/aer/ibmq_something
    'device':'ibm_nairobi' #to be simulated if using simulator
    }
algorithm = Algorithm(options)
optimizer = Optimizer('spsa', disp=False, max_iter=100)
vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(ising=True), optimizer, algorithm)
"""
import time
start_time = time.time()



#ground state
optimized_parameters = vqe.optimize_parameters(vqe.vqe_expval)
eigenvalue = vqe.vqe_expval(optimized_parameters)
eigenstate = vqe.get_eigenstate(optimized_parameters)
print("result: ", eigenvalue, " expected: ", (-B/2-J/4)*orbitals, eigenstate)

print("--- %s seconds ---" % (time.time() - start_time))

"""


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