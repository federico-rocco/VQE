# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 12:18:37 2022

@author: cosmo
"""

from qiskit import QuantumCircuit, transpile
from qiskit.tools.monitor import job_monitor
from qiskit.providers.aer.noise import NoiseModel

import qiskit as qk
import numpy as np


class Algorithm:
    #handles execution
    
    def __init__(self, options):
        self.shots = 1024 if options.get('shots') == None\
                            else options.get('shots')
        self.seed = np.random.randint(0,100000) if options.get('seed') == None\
                            else options.get('seed')        
        self.ibmq = options.get('ibmq')        

        
        
        if self.ibmq == True:
            provider = qk.IBMQ.load_account()
            if options.get('backend')=='least busy':
                from qiskit.providers.ibmq import least_busy
                self.backend = least_busy(provider.backends(filters=lambda x: x.configuration().n_qubits >= 3 and not x.configuration().simulator and x.status().operational==True))
            else:
                self.backend = provider.get_backend(options.get('backend'))
            self.monitor = job_monitor

            
        else:
            if options.get('backend') == None:
                options['backend'] = 'qasm_simulator' 
            self.backend = qk.Aer.get_backend(options['backend'])            
            

            
            self.noise_model, self.coupling_map, self.basis_gates = None, None, None
            if options.get('device') != None: #simulate a real device
                device = options.get('device')
                provider = qk.IBMQ.load_account()
                backend = provider.get_backend(device)                
                self.noise_model = NoiseModel.from_backend(backend)
                self.coupling_map = backend.configuration().coupling_map
                self.basis_gates = self.noise_model.basis_gates

    
    def measure(self, qc, qb, cb, n=0):
        
        #qc = transpile(qc, qk.Aer.get_backend('aer_simulator_statevector'), optimization_level=3)
   
        qc.measure(qb,cb)
        
        if self.ibmq:
            job = qk.execute(qc, 
                            backend = self.backend, 
                            shots=self.shots)
            self.monitor(job)
            
        else:
            job = qk.execute(qc, 
                            backend = self.backend, 
                            shots=self.shots,
                            seed_transpiler=self.seed,
                            seed_simulator=self.seed,
                            noise_model=self.noise_model,
                            coupling_map=self.coupling_map,
                            basis_gates=self.basis_gates)
            """
                            optimization_level=self.optimization_level,
                            ,
                            
                            initial_layout=self.layout)"""

        return job.result().get_counts(qc)

        
        
        
        
        

        
        
        
        
        
        
        
        
        
        
        
        