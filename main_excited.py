# -*- coding: utf-8 -*-
"""
Created on Wed Jul 13 18:49:15 2022

@author: cosmo
"""

from hamiltonian import Hamiltonian
from ansatz import Ansatz
from solver import Eigensolver, State
from qiskit import Aer, IBMQ
from algorithm import Algorithm
from optimizer import Optimizer
from quarkoniumpaper import *
import matplotlib.pyplot as plt
import matplotlib as mpl


fermions = 1
orbitals = 3
model = Quarkonium()
coeffs = model.coeff(orbitals)
hamiltonian = Hamiltonian(fermions,orbitals,coeffs)
ansatz = Ansatz('quarkonium', n_fermions=fermions, n_qubits=orbitals, mp2=False)

options = {
    'shots':1024*100,
    'ibmq':False,
    'backend':'qasm_simulator', #qasm/aer/ibmq_something
    #'device':'ibm_nairobi' #to be simulated if using simulator
    }
algorithm = Algorithm(options)
optimizer = Optimizer('cobyla', disp=True, max_iter=100)



x = []
y = []
en = []
for delta in range(1000,1001,20):
    x.append(delta)
    wrong = 0
    wrong_eig = []
    for i in range(1):       
        vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(), optimizer, algorithm)
        vqe.factor = delta 
        #ground state
        optimized_parameters = vqe.optimize_parameters(vqe.analytic_expval)
        eigenvalue = vqe.analytic_expval(optimized_parameters)
        eigenstate = vqe.get_eigenstate(optimized_parameters)
        vqe.save_state(State(eigenvalue, eigenstate, optimized_parameters))
        en.append(eigenvalue)
        print(eigenvalue)
        print(eigenstate)
        
        #first excited
        optimized_parameters1 = vqe.optimize_parameters(vqe.expval_excited_state)
        eigenvalue1 = vqe.analytic_expval(optimized_parameters1)
        eigenstate1 = vqe.get_eigenstate(optimized_parameters1)
        vqe.save_state(State(eigenvalue1, eigenstate1, optimized_parameters1))
        en.append(eigenvalue1)
        print(eigenvalue1)
        print(eigenstate1)
        
        
        #second excited
        optimized_parameters2 = vqe.optimize_parameters(vqe.expval_excited_state)
        eigenvalue2 = vqe.analytic_expval(optimized_parameters2)
        eigenstate2 = vqe.get_eigenstate(optimized_parameters2)
        vqe.save_state(State(eigenvalue2, eigenstate2, optimized_parameters2))
        en.append(eigenvalue2)
        print(eigenvalue2)
        print(eigenstate2)
        

        if not eigenvalue2 > eigenvalue1 > eigenvalue:
            wrong += 1
            wrong_eig.append(eigenvalue)
            wrong_eig.append(eigenvalue1)
            wrong_eig.append(eigenvalue2)
            wrong_eig.append("---")

     
    y.append(wrong*5)       
    print("With delta/", delta, wrong*5, "% error.")
    print("========================")
    
  

#en = [3123.8625522666325-2*1317, 3849.9106274974165-2*1317, 4777.029643873779-2*1317]
en.sort()
orbital=[r'1S',r'2S',r'3S']
x = [1] * len(en)

fig, ax = plt.subplots()
ax.scatter(x, en, s=9000, marker="_", linewidth=2, zorder=3)
ax.grid(axis='y')


for xi,yi,tx in zip(x,en,orbital):
    ax.annotate(tx, xy=(xi,yi), xytext=(7,-3), size=10,
                ha="center",va='top', textcoords="offset points")
    
ax.margins(0.1)
plt.ylabel('Energy [MeV]')
plt.title('Charmonium Energy levels')
ax.set_xticks([])
ax.yaxis.set_minor_locator(mpl.ticker.MaxNLocator(50))
#plt.plot(x,y)
#plt.show