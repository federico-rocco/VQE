# -*- coding: utf-8 -*-
"""
Created on Fri Sep 23 16:12:42 2022

@author: cosmo
"""

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
import os

coeffs = coeff()
hamiltonian = Hamiltonian(fermions,orbitals,coeffs)
ansatz = ansatz('quarkonium', n_fermions=fermions, n_qubits=orbitals)


options = {
    'seed':1,
    'shots':1024,
    'ibmq':False,
    'backend':'qasm_simulator', #qasm/aer/ibmq_something
    #'device':'ibmq_athens' #to be simulated if using simulator
    }
algorithm = Algorithm(options)
optimizer = Optimizer('spsa', max_iter=100)
vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(), optimizer, algorithm)
optimized_parameters = [3.31,0.95]
energy = vqe.vqe_expval()
print(energy)