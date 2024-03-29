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
from qiskit.providers.aer import noise

import pickle
import os

class Algorithm:
    #handles execution
    
    def __init__(self, options, n_qubits=3):
        self.shots = 1024 if options.get('shots') == None\
                            else options.get('shots')
        self.seed = None if options.get('seed') == None\
                            else options.get('seed')        
        self.ibmq = options.get('ibmq') 
        self.noise_model, self.coupling_map, self.basis_gates = None, None, None
        self.n_qubits = n_qubits

        
        
        if self.ibmq == True:
            provider = qk.IBMQ.load_account()
            if options.get('backend')=='least busy':
                from qiskit.providers.ibmq import least_busy
                self.backend = least_busy(provider.backends(filters=lambda x: x.configuration().n_qubits >= self.n_qubits and not x.configuration().simulator and x.status().operational==True))
            else:
                self.backend = provider.get_backend(options.get('backend'))
            self.monitor = job_monitor

            
        else:
            if options.get('backend') == None:
                options['backend'] = 'qasm_simulator' 
            self.backend = qk.Aer.get_backend(options['backend'])

            
            
            if options.get('device') != None: #simulate a real device
                provider = qk.IBMQ.load_account()
                name_device = options.get('device')
                
                try:
                    __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
                    self.noise_model = NoiseModel.from_dict(pickle.load(open(__location__+'/noise_models/'+name_device+'.pkl','rb')))
                    self.basis_gates = self.noise_model.basis_gates
                    
                except FileNotFoundError:
                    device = provider.get_backend(name_device)                   
                    self.noise_model = NoiseModel.from_backend(device)
                    self.backend_properties = device.properties()
                    self.coupling_map = device.configuration().coupling_map
                    self.basis_gates = self.noise_model.basis_gates
                    self.transpiler = device


    def measure(self, qc, qb, cb):
        qc.measure(qb, cb)
            
        if self.noise_model is not None:
            qc = transpile(qc,   
                           backend=self.backend,
                           #backend_properties=self.backend_properties,
                           #coupling_map=self.coupling_map,
                           basis_gates=self.basis_gates)
            #print(qc)
            
        
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
                            noise_model=self.noise_model)
                            #coupling_map=self.coupling_map,
                            #basis_gates=self.basis_gates)
            """
                            optimization_level=self.optimization_level,
                            ,
                            
                            initial_layout=self.layout)"""

        return job.result().get_counts(qc)

        
        
        
        
        

        
        
        
        
        
        
        
        
        
        
        
        