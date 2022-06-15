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
from qiskit import Aer, QuantumCircuit
from qiskit.algorithms import VQE

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


hamiltonian = (1/2*(21/4*omega + V(0,0) + V(1,1) + V(2,2)) * I ^ I ^ I) + \
        (-1/2*(3/4*omega + V(0,0)) * I ^ I ^ Z) + \
        (-1/2*(7/4*omega + V(1,1)) * I ^ Z ^ I) + \
        (-1/2*(11/4*omega + V(2,2)) * Z ^ I ^ I) + \
        (1/4*(-np.sqrt(3/2)*omega + 2*V(0,1)) * I ^ X ^ X) + \
        (1/4*(-np.sqrt(5)*omega + 2*V(1,2)) * X ^ X ^ I) + \
        (1/4*(-np.sqrt(3/2)*omega + 2*V(0,1)) * I ^ Y ^ Y) + \
        (1/4*(-np.sqrt(5)*omega + 2*V(1,2)) * Y ^ Y ^ I) + \
        (1/2*V(0,2) * X ^ Z ^ X) + \
        (1/2*V(0,2) * Y ^ Z ^ Y) 
        
num_particles = [0,1]
num_spin_orbitals = 3


qubit_converter = QubitConverter(mapper=JordanWignerMapper())


#qubit_op = qubit_converter.convert(hamiltonian,num_particles=num_particles)

qc = QuantumCircuit(3,3)
"""
qc.ry(beta, 1)
qc.ry(2*alpha, 2)
qc.cx(2,0)
qc.cx(0,1)
qc.x(2)
qc.ry(-beta,1)
qc.cx(0,1)
qc.cx(1,0)
"""
initial_state = qc

reps=1
var_form = UCC(
    excitations="sd",
    num_particles=num_particles,
    num_spin_orbitals=num_spin_orbitals,
    initial_state=initial_state,
    qubit_converter=qubit_converter,
    reps=reps,
)
optimizer=COBYLA(maxiter=1000)
vqe = VQE(ansatz=var_form, optimizer=optimizer, quantum_instance=Aer.get_backend("aer_simulator_statevector"))
vqe_result = vqe.compute_minimum_eigenvalue(hamiltonian)
print(vqe_result)


    
    
    
    
    
    
    
    
    
    
    
    