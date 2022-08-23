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
from optimizer import Minimizer
import numpy as np 

def num_integrate_gs(B):
    """
    numerically integrate exact band to get gs energy of TIM
    this should give -E_0/(N*J) by Pfeufy
    Here set J=1 (units of energy)
    """
    # lamba_ratio (setting J=1): compare thesis
    ll = 1/(2*B)
    
    # set energy
    gs_energy = 0
    
    # numerical integration
    step_size = 0.0001
    k_values = np.arange(0, np.pi, step_size)
    integration_values = [step_size*np.sqrt(1 + ll**2 + 2*ll*np.cos(kk)) for kk in k_values]
    integral = np.sum(integration_values)
    gs_energy = 1*integral/(4*np.pi*ll)
    
    return gs_energy

orbitals = 3
fermions = 0
B = 0.1
ising_coeff = {'B':B,
               'J':1
               }
hamiltonian = Hamiltonian(n_qubits=orbitals, coeff=ising_coeff)
ansatz = ansatz('ising', n_qubits=orbitals)

options = {
    'shots':1024,
    'ibmq':False,
    'seed':10,
    'backend':'aer_simulator_statevector', #qasm/aer/ibmq_something
    #'device':'ibm_nairobi' #to be simulated if using simulator
    }
algorithm = Algorithm(options)
optimizer = Minimizer('spsa', disp=False, max_iter=100)
vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(ising=True), optimizer, algorithm)

import time
start_time = time.time()



#ground state
optimized_parameters = vqe.optimize_parameters(vqe.vqe_expval)
eigenvalue = vqe.vqe_expval(optimized_parameters)
print("result: ", eigenvalue, num_integrate_gs(B))

print("--- %s seconds ---" % (time.time() - start_time))