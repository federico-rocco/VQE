# -*- coding: utf-8 -*-
"""
Created on Thu May 12 10:59:48 2022

@author: cosmo
"""

from qiskit import IBMQ
from qiskit.providers.ibmq import least_busy
from ferm_op import Hamiltonian


provider = IBMQ.load_account();

state = Hamiltonian(1,3)
hamiltonian = state.buildH()
ansatz = state.build_ansatz

num_qubits = 3
backend = least_busy(provider.backends(filters=lambda x: x.configuration().n_qubits >= num_qubits and not x.configuration().simulator and x.status().operational==True))
from qiskit.tools.monitor import job_monitor
job_monitor(job)
vqe = state.vqe(backend)


# This is the target energy
h2_energy = 492.6

result = vqe.compute_minimum_eigenvalue(hamiltonian)
print('Result:', result.optimal_value, 'Reference:', h2_energy)
print(result)