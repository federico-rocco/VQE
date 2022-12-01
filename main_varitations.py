# -*- coding: utf-8 -*-
"""
Created on Mon Oct 31 09:58:23 2022

@author: cosmo
"""

from hamiltonian import Hamiltonian
from ansatz import Ansatz
from solver import Eigensolver, State
from qiskit import Aer, IBMQ
from algorithm import Algorithm
from optimizer import Optimizer
from quarkoniumpaper import Quarkonium
import matplotlib as mpt

fermions = 1
orbitals = 3

model = Quarkonium()
coeffs = model.coeff(orbitals)
hamiltonian = Hamiltonian(fermions,orbitals,coeffs)



options = {
    'shots':1024*1000,
    'ibmq':False,
    'seed':1,
    'backend':'qasm_simulator', #qasm/aer/ibmq_something
    #'device':'ibm_nairobi' #to be simulated if using simulator
    }
algorithm = Algorithm(options)
ansatz = Ansatz('quarkonium', n_fermions=fermions, n_qubits=orbitals, mp2=False)

optimizers = ['spsa']
for opt in optimizers:
    optimizer = Optimizer(opt, disp=True, max_iter=1000000)
    
    vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(), optimizer, algorithm)

    import time
    start_time = time.time()
    
    
    
    #ground state
    optimized_parameters = vqe.optimizer.opt.minimize(vqe.analytic_expval,[0,0])
    
    #eigenvalue = vqe.expval(optimized_parameters)
    print(opt)
    print(optimized_parameters)
    print("--- %s seconds ---" % (time.time() - start_time))

    #print(vqe.get_eigenstate(optimized_parameters.x))
    print('-----')
    
#mpt.pyplot.hist(eigs)
