# -*- coding: utf-8 -*-
"""
Created on Wed Jul 20 19:22:37 2022

@author: cosmo
"""

from hamiltonian import Hamiltonian
from ansatz import Ansatz
from solver import Eigensolver, State
from qiskit import Aer, IBMQ
from algorithm import Algorithm
from optimizer import Optimizer
from quarkonium import Quarkonium
import matplotlib as mpt

fermions = 1
orbitals = 3
model = Quarkonium('charmonium')
coeffs = model.coeff(orbitals)
hamiltonian = Hamiltonian(fermions,orbitals,coeffs)
ansatz = Ansatz('quarkonium', n_fermions=fermions, n_qubits=orbitals, mp2=False)

options = {
    'shots':1024*100,
    'ibmq':False,
    'seed':1,
    'backend':'qasm_simulator', #qasm/aer/ibmq_something
    #'device':'ibm_nairobi' #to be simulated if using simulator
    }
algorithm = Algorithm(options)
optimizer = Optimizer('spsa', disp=False, max_iter=200)

vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(), optimizer, algorithm)

import time
start_time = time.time()

eigs = []
for i in range(1):
    #ground state
    optimized_parameters = vqe.optimize_parameters(vqe.analytic_expval)
    
    eigenvalue = vqe.expval(optimized_parameters)
    print("spin averaged:", eigenvalue + 2*model.mq)
    print("singlet:", eigenvalue + model.VSS(0), "triplet:", eigenvalue + model.VSS(1))
    print("ηc:", eigenvalue + model.VSS(0) + 2*model.mq, "J/Ψ:", eigenvalue + model.VSS(1) + 2*model.mq)
    eigs.append(eigenvalue)
    

#mpt.pyplot.hist(eigs)
print("--- %s seconds ---" % (time.time() - start_time))

print(model.VSS(1)-model.VSS(0))