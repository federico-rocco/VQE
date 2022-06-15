# -*- coding: utf-8 -*-
"""
Created on Thu May 12 11:45:00 2022

@author: cosmo
"""

from qiskit_nature.settings import settings
from qiskit_nature.drivers import UnitsType
#from qiskit_nature.drivers.second_quantization import JordanWignerMapper, ParityMapper, FermionicMapper
from qiskit_nature.mappers.second_quantization import JordanWignerMapper, ParityMapper, FermionicMapper
from qiskit_nature.problems.second_quantization.electronic import ElectronicStructureProblem
from qiskit_nature.converters.second_quantization import QubitConverter
from qiskit_nature.operators.second_quantization import FermionicOp
from qiskit_nature.circuit.library import HartreeFock, UCC
from qiskit.algorithms import VQE
from qiskit.algorithms.optimizers import COBYLA, SPSA, SLSQP, CG
from qiskit import Aer, transpile
from qiskit.utils import QuantumInstance
from qiskit.opflow import I, X, Z, Y, StateFn
from qiskit.opflow.gradients import Gradient
from qiskit.circuit import QuantumCircuit, Parameter
from qiskit.circuit.library import TwoLocal
from qiskit.quantum_info import state_fidelity, Statevector

import random
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import mpmath

k=0.4063
sigma=441.6**2
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


def coeff(m,n):
    return T(m,n) + V(m,n)
        
class Hamiltonian:
    
    def __init__(self, n_qubit):    
        self.n_qubit = n_qubit
        self.wavefunction = []
        self.hamiltonian = []
        self.alpha = Parameter('alpha')
        self.beta = Parameter('beta')
    
    
    def build_a(self, i, dag=True):
        label = ""
        for item in range(self.n_qubit):        
            if item == i:
                if dag==True:
                    label += "+"
                else:
                    label += "-"
            else:
                label += "I" 
        return label
    
    
    def buildH(self):
        ferm_op = 0
        for m, n in [[_m, _n] for _m in range(self.n_qubit) for _n in range(self.n_qubit)]:
            ferm_op += FermionicOp((self.build_a(m,True),coeff(m,n)), display_format="dense") @ FermionicOp(self.build_a(n,False), display_format="dense")            
        hamiltonian = JordanWignerMapper().map(ferm_op)
        self.hamiltonian = hamiltonian
        return hamiltonian
    
    def build_ansatz(self):
        hartree = HartreeFock(3, (1, 0), qubit_converter=QubitConverter(mapper=JordanWignerMapper()))
        #wavefunction = hartree
        #print(hartree)
        hartree = QuantumCircuit(3)
        wavefunction = hartree
        #wavefunction.x(0)
        wavefunction.ry(self.beta, 1)
        wavefunction.ry(2*self.alpha, 2)
        wavefunction.cx(2, 0)
        wavefunction.cx(0, 1)
        wavefunction.x(2)
        wavefunction.ry(-1*self.beta, 1)
        wavefunction.cx(0, 1)
        wavefunction.cx(1, 0)
        self.wavefunction = wavefunction
        print(wavefunction)
        return wavefunction


    def vqe(self, backend):
        hamiltonian = self.buildH()
        wavefunction = self.build_ansatz()
        op = ~StateFn(hamiltonian)@StateFn(wavefunction)
        grad = Gradient(grad_method='lin_comb').convert(operator=op, params=[self.alpha,self.beta])
        qi_sv = QuantumInstance(backend,
                        shots=1000,
                        seed_simulator=2,
                        seed_transpiler=2)
        optimizer = SPSA(maxiter=100)        
        vqe = VQE(wavefunction, optimizer=optimizer, gradient=grad, quantum_instance=qi_sv)
        return vqe
"""
state = Hamiltonian(3)
h = state.buildH()
a = state.build_ansatz()
op = ~StateFn(h)@StateFn(a)
print(op)
grad = Gradient(grad_method='lin_comb').convert(operator=op, params=[state.alpha,state.beta])
"""




