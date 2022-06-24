# -*- coding: utf-8 -*-
"""
Created on Sat Apr 30 19:08:47 2022

@author: cosmo
"""


from qiskit.opflow import X, Y, Z, I
from qiskit_nature.settings import settings
from qiskit_nature.drivers import UnitsType
from qiskit_nature.drivers.second_quantization import PySCFDriver
from qiskit_nature.problems.second_quantization.electronic import ElectronicStructureProblem
from qiskit_nature.mappers.second_quantization import ParityMapper
from qiskit_nature.converters.second_quantization import QubitConverter
from qiskit_nature.mappers.second_quantization import JordanWignerMapper
from qiskit_nature.circuit.library import HartreeFock, UCC
from qiskit.algorithms.optimizers import COBYLA, SPSA, SLSQP
from qiskit import Aer, QuantumCircuit, QuantumRegister
from qiskit.algorithms import VQE
from qiskit.circuit import Parameter

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import mpmath


theta = Parameter('theta')
n_fermi = 1
n_qubits = 3
qb = QuantumRegister(n_qubits)
qc = QuantumCircuit(qb)
for i in range(n_fermi):
            for a in range(n_fermi,n_qubits):
                
                # Pauli-z
                #for k in range(i+1,a):
                #    #qc.z(qb[k])
                #    qc.cx(qb[k],qb[a]) 

                # Pauli strings
                qc.h(qb[i])
                qc.rx(np.pi/2,qb[i])
                qc.h(qb[a])

                qc.cx(qb[i],qb[a])

                #qc.rz(theta[i,a-n_fermi]/2,qb[a])
                qc.rz(theta/2,qb[a])

                qc.cx(qb[i],qb[a])

                ## Pauli-z
                #for k in range(i+1,a):
                #    #qc.z(qb[k])
                #    qc.cx(qb[k],qb[a]) 

                qc.rx(-np.pi/2,qb[i])
                qc.rx(np.pi/2,qb[a])

                # Pauli-z
                #for k in range(i+1,a):
                #    #qc.z(qb[k])
                #    qc.cx(qb[k],qb[a]) 

                qc.cx(qb[i],qb[a])

                #qc.rz(-theta[i,a-n_fermi]/2,qb[a])
                qc.rz(-theta/2,qb[a])

                qc.cx(qb[i],qb[a])

                #for k in range(i+1,a):
                #    #qc.z(qb[k])
                #    qc.cx(qb[k],qb[a])

                qc.h(qb[i])
                qc.rx(-np.pi/2,qb[a])
                qc.h(qb[a])
                
qc.draw('mpl')


    
    
    
    
    
    
    
    
    
    
    
    