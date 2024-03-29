# -*- coding: utf-8 -*-
"""
Created on Thu May 12 10:59:48 2022

@author: cosmo
"""

from hamiltonian import Hamiltonian
from ansatz import Ansatz
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


fermions = 1
orbitals = 3
model = Quarkonium('charmonium')
coeffs = model.coeff(orbitals)
hamiltonian = Hamiltonian(fermions,orbitals,coeffs)
ansatz = Ansatz('quarkonium', n_fermions=fermions, n_qubits=orbitals)


options = {
    #'seed':1,
    'shots':20000,
    'ibmq':False,
    'backend':'qasm_simulator', #qasm/aer/ibmq_something
    'device':'ibm_nairobi' #to be simulated if using simulator
    }
algorithm = Algorithm(options)
optimizer = Optimizer('cobyla', max_iter=1000)
vqe = Eigensolver(fermions, orbitals, ansatz, hamiltonian(), optimizer, algorithm)
optimized_parameters = [1.43012171, 6.47848588]#vqe.optimize_parameters(vqe.expval)

print('found parameters')

executions=100

x = []
y = []
y_std= []
weights = []
for lamda in range(1,7,1):
    x.append(lamda)
    print("lamda = ", lamda)
    vqe.set_folding(lamda)
    energies = []
    for i in range(executions):
        
        energy = vqe.expval(optimized_parameters)
        energies.append(energy)
    mean_energy = np.mean(energies)
    std = np.std(energies)
    y.append(mean_energy)
    y_std.append(std)
    #weights.append(1/std)
    print("Optimized energy: ", mean_energy)
    
    print("------------------------")

x = np.array(x)
y = np.array(y)
b, m = np.polynomial.polynomial.polyfit(x, y, 1, rcond=None, full=False)
#plt.errorbar(x,y,yerr=y_std)
#plt.xlim([-1,11])
#plt.ylim([450,700])

x2 = np.insert(x, 0, 0)

y2 = np.insert(y, 0, b)
y_std = np.array(y_std)
y_std2 = np.insert(y_std, 0, np.mean(y_std))

plt.errorbar(x,y,yerr=y_std, ecolor='tab:red', capsize=3, fmt="r--o",linewidth=1, mfc='none')
plt.plot(x2, y2, linestyle='dashed', color='red')
plt.plot([0],[b],'r*')
plt.fill_between(x2, y2-y_std2, y2+y_std2, edgecolor='pink', facecolor='moccasin')
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