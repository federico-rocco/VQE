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
from quarkoniumpaper import Quarkonium
import matplotlib as mpt
from scipy.optimize import root
import numpy as np

def func(theta):
    return [vqe.analytic_expval(theta)-492.6, 0]

fermions = 1
orbitals = 3
flavours = ['charmonium']
for flavour in flavours:
    model = Quarkonium(flavour)
    coeffs = model.coeff(orbitals)
    hamiltonian = Hamiltonian(fermions,orbitals,coeffs)
    ansatz = Ansatz('UCCSD', n_fermions=fermions, n_qubits=orbitals, mp2=False)
    
    options = {
        'shots':1024*100,
        'ibmq':False,
        'seed':1,
        'backend':'qasm_simulator', #qasm/aer/ibmq_something
        #'device':'ibm_nairobi' #to be simulated if using simulator
        }
    algorithm = Algorithm(options)
    optimizer = Optimizer('cobyla', disp=True, max_iter=2000)
    
    vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(), optimizer, algorithm)
    
    import time
    start_time = time.time()
    #initial_guess=np.array([0,0])
    #optimized_parameters = root(func, initial_guess).x

    #ground state
    optimized_parameters = vqe.optimize_parameters(vqe.analytic_expval)
    
    eigenvalue = vqe.expval(optimized_parameters)
    print(eigenvalue, vqe.get_eigenstate(optimized_parameters))
    #print("eigenstate: ", vqe.get_eigenstate(optimized_parameters))
    """
    print(flavour)
    print("spin averaged:", eigenvalue + 2*model.mq)
    #print("singlet:", eigenvalue + model.VSS(0), "triplet:", eigenvalue + model.VSS(1))
    print("singlet:", eigenvalue + model.VSS(0) + 2*model.mq, "triplet:", eigenvalue + model.VSS(1) + 2*model.mq)
    
    """
    print('--------')
    
    #mpt.pyplot.hist(eigs)
    print("--- %s seconds ---" % (time.time() - start_time))
