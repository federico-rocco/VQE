# -*- coding: utf-8 -*-
"""
Created on Mon May 30 14:38:20 2022

@author: cosmo
"""

from qiskit import Aer
from qiskit.utils import QuantumInstance

# Algorithm Imports
from qiskit.algorithms import VQE, NumPyMinimumEigensolver
from qiskit.algorithms.optimizers import CG

from qiskit.opflow import I, X, Z, Y, StateFn
from qiskit.circuit import QuantumCircuit, ParameterVector
from scipy.optimize import minimize
from qiskit.opflow.gradients import Gradient
from qiskit_nature.operators.second_quantization import FermionicOp
from qiskit_nature.mappers.second_quantization import JordanWignerMapper, ParityMapper, FermionicMapper

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import mpmath


k=0.4063
sigma=441.6
mu=637.5
omega=562.9
b=1/np.sqrt(mu*omega)

def kron_delta(a,b):
    if a == b :
        return 1
    else:
        return 0
    
def T(m,n):
    return (omega/2*((2*n+3/2)*kron_delta(m,n)-
           np.sqrt(n*(n+1/2))*kron_delta(m+1,n)-
           np.sqrt((n+1)*(n+3/2))*kron_delta(m-1,n)))

def r(m,n):
    return ((-1)**(m+n)*4*b/(np.pi*(1-4*n**2))*np.sqrt(sp.special.gamma(m+3/2)*
           sp.special.gamma(n+3/2)/(sp.special.factorial(m)*sp.special.factorial(n)))*
           float(mpmath.hyp2f1(2,-m,3/2-n,1)))

def r_inv(m,n):
    return ((-1)**(m+n)*4/(b*np.pi*(1+2*n))*np.sqrt(sp.special.gamma(m+3/2)*
           sp.special.gamma(n+3/2)/(sp.special.factorial(m)*sp.special.factorial(n)))*
           float(mpmath.hyp3f2(1/2,1,-m,3/2,1/2-n,1)))

def V(m,n):
    return -k*r_inv(m,n)+sigma*r(m,n)

# Instantiate the system Hamiltonian
h2_hamiltonian = (1/2*(21/4*omega + V(0,0) + V(1,1) + V(2,2)) * I ^ I ^ I) + \
        (-1/2*(3/4*omega + V(0,0)) * I ^ I ^ Z) + \
        (-1/2*(7/4*omega + V(1,1)) * I ^ Z ^ I) + \
        (-1/2*(11/4*omega + V(2,2)) * Z ^ I ^ I) + \
        (1/4*(-np.sqrt(3/2)*omega + 2*V(0,1)) * I ^ X ^ X) + \
        (1/4*(-np.sqrt(5)*omega + 2*V(1,2)) * X ^ X ^ I) + \
        (1/4*(-np.sqrt(3/2)*omega + 2*V(0,1)) * I ^ Y ^ Y) + \
        (1/4*(-np.sqrt(5)*omega + 2*V(1,2)) * Y ^ Y ^ I) + \
        (1/2*V(0,2) * X ^ Z ^ X) + \
        (1/2*V(0,2) * Y ^ Z ^ Y) 
        


# This is the target energy
h2_energy = 492.6

# Define the Ansatz
wavefunction = QuantumCircuit(3)
params = ParameterVector('theta', length=2)
it = iter(params)
wavefunction.ry(params[1], 1)
wavefunction.ry(2*params[0], 2)
wavefunction.cx(2, 0)
wavefunction.cx(0, 1)
wavefunction.x(2)
wavefunction.ry(-1*params[1], 1)
wavefunction.cx(0, 1)
wavefunction.cx(1, 0)

# Define the expectation value corresponding to the energy
op = ~StateFn(h2_hamiltonian) @ StateFn(wavefunction)

grad = Gradient(grad_method='lin_comb')

qi_sv = QuantumInstance(Aer.get_backend('aer_simulator_statevector'),
                        shots=1,
                        seed_simulator=2,
                        seed_transpiler=2)

#Conjugate Gradient algorithm
optimizer = CG(maxiter=50)

# Gradient callable
vqe = VQE(wavefunction, optimizer=optimizer, gradient=grad, quantum_instance=qi_sv)

result = vqe.compute_minimum_eigenvalue(h2_hamiltonian)
print('Result:', result.optimal_value, 'Reference:', h2_energy)