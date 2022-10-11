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

import time
start_time = time.time()


coeffs = coeff('charmonium')
hamiltonian = Hamiltonian(fermions,orbitals,coeffs)
ansatz = ansatz('quarkonium', n_fermions=fermions, n_qubits=orbitals)


options = {
    'seed':1,
    'shots':1024*10,
    'ibmq':False,
    'backend':'qasm_simulator', #qasm/aer/ibmq_something
    'device':'ibm_nairobi' #to be simulated if using simulator
    }
algorithm = Algorithm(options)
optimizer = Optimizer('spsa', max_iter=200)
vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(), optimizer, algorithm)
optimized_parameters = [1.65778537, 3.27087804]#vqe.optimize_parameters(vqe.expval)



executions=1

x = []
y = []
y_std= []
weights = []
for lamda in range(1,6,1):
    x.append(lamda)
    print("lamda = ", lamda)
    vqe.set_folding(lamda)
    energies = []
    for i in range(executions):
        
        energy = vqe.expval(optimized_parameters)
        energies.append(energy)
    mean_energy = np.mean(energies)
    #std = np.std(energies)
    y.append(mean_energy)
    #y_std.append(std)
    #weights.append(1/std)
    print("Optimized energy: ", mean_energy)
    
    print("------------------------")

b, m = np.polynomial.polynomial.polyfit(x, y, 1, rcond=None, full=False)
#plt.errorbar(x,y,yerr=y_std)
#plt.xlim([-1,11])
#plt.ylim([450,700])
x = np.array(x)
plt.plot(x, m*x + b)
plt.scatter(x, y)
plt.show()


print("--- %s seconds ---" % (time.time() - start_time))
"""
script_dir = os.path.dirname(__file__)
output_dir = os.path.join(script_dir, "GraphicsOutput/Noise/")
if not os.path.isdir(output_dir):
    os.makedirs(output_dir)
    
#plot
fig,ax = plt.subplots()
plt.rc('font', family='monospace')
plt.text(0.5, 1.07, "Zero noise extrapolation", horizontalalignment='center', fontsize=12, transform = ax.transAxes)
ax.plot(x, y, color="blue", linestyle="-", linewidth=1, label = 'P Newton')
ax.plot(x, m*x+b, color="black", linestyle="-", linewidth=2,  label = 'P tov')
ax.set_xlabel('λ',fontsize=14)
ax.set_ylabel('E [MeV]', fontsize=14)
ax.minorticks_on()

ax2 = ax.twinx()
ax2.minorticks_on()
ax2.plot(r_newton, m_newton,color="blue", linestyle=":", label = 'm Newton')
ax2.plot(r_tov, m_tov, color="black", linestyle="-.", label = 'm tov')
ax2.plot(R_newton, M_newton, marker = 'o', linestyle="", color='green', label='NS Newton mass')
ax2.plot(R_tov, M_tov, marker = 'o', color='red', linestyle="", label='NS tov mass')
ax2.set_ylabel(r"m [$M_{\odot}$]",fontsize=14)

fig.legend(loc='center right', bbox_to_anchor=(1,1), bbox_transform=ax.transAxes)
plt.rcParams["savefig.bbox"] = "tight"

fig.savefig(output_dir+'zero_noise_extrapolation.pdf', format='pdf', dpi=1000)
"""