# -*- coding: utf-8 -*-
"""
Created on Thu May 12 10:59:48 2022

@author: cosmo
"""

from qiskit import IBMQ, QuantumCircuit, execute
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
from qiskit import Aer
from qiskit.algorithms import VQE


import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import mpmath
from ferm_op import Hamiltonian

options = {
	'backend_name': 'ibmq_qasm_simulator'
}


state = Hamiltonian(3)
hamiltonian = state.buildH()
ansatz = state.build_ansatz

backend = Aer.get_backend('aer_simulator_statevector')

vqe = state.vqe(backend)


# This is the target energy
h2_energy = 492.6

result = vqe.compute_minimum_eigenvalue(hamiltonian)
#print('Result:', result.optimal_value, 'Reference:', h2_energy)
#print(result)