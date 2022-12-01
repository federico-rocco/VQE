# -*- coding: utf-8 -*-
"""
Created on Wed Oct 19 16:40:12 2022

@author: cosmo
"""

from hamiltonian import Hamiltonian
from ansatz import Ansatz
from solver import Eigensolver, State
from qiskit import Aer, IBMQ
from algorithm import Algorithm
from optimizer import Optimizer
from quarkonium2 import Quarkonium
import matplotlib as mpt
import matplotlib.pyplot as plt

fermions = 1
orbitals = 3
import time
start_time = time.time()
x=[]
y=[]
for w in range (160,900,20):
    model = Quarkonium('bc', w)
    coeffs = model.coeff(orbitals)
    hamiltonian = Hamiltonian(fermions,orbitals,coeffs)
    ansatz = Ansatz('quarkonium', n_fermions=fermions, n_qubits=orbitals, mp2=False)
    
    options = {
        'shots':1024*10,
        'ibmq':False,
        'seed':1,
        'backend':'qasm_simulator', #qasm/aer/ibmq_something
        #'device':'ibm_nairobi' #to be simulated if using simulator
        }
    algorithm = Algorithm(options)
    optimizer = Optimizer('slsqp', disp=False, max_iter=1000)
    
    vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(), optimizer, algorithm)

    optimized_parameters = vqe.optimize_parameters(vqe.analytic_expval)
    
    eigenvalue = vqe.expval(optimized_parameters)
    print("omega: ", w, "energy:", eigenvalue)
    if True:
        x.append(w)
        y.append(eigenvalue)


fig, ax = plt.subplots()
ax.scatter(x, y, marker=".", linewidth=2, zorder=3)
#ax.grid(axis='y')

    
#ax.margins(0.1)
plt.ylabel('Energy [MeV]')
plt.xlabel('Frequency [MeV]')

    
#mpt.pyplot.scatter(x,y)
#mpt.pyplot.hist(eigs)
print("--- %s seconds ---" % (time.time() - start_time))