# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 15:33:26 2022

@author: cosmo
"""

from hamiltonian import Hamiltonian
from ansatz import ansatz
from solver import Eigensolver, State
from qiskit import Aer, IBMQ, QuantumCircuit
from algorithm import Algorithm
from optimizer import Optimizer
from quarkonium import *
import matplotlib as mpt

coeffs = coeff()
hamiltonian = Hamiltonian(fermions,orbitals,coeffs)
ansatz = ansatz('', n_fermions=fermions, n_qubits=orbitals)

options = {
    'shots':1024,
    'ibmq':True,
    'seed':10,
    'backend':'ibm_oslo', #qasm/aer/ibmq_something
    #'device':'ibm_nairobi' #to be simulated if using simulator
    }
algorithm = Algorithm(options)
optimizer = Optimizer('', disp=False, max_iter=10)

vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(), optimizer, algorithm)
qc = QuantumCircuit(3)
qc.h(0)
qc.x(1)
qc.y(1)
qc.z(1)
vqe.set_folding(2)
qc = vqe.fold(qc)
print(qc)
