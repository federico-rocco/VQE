# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 12:18:37 2022

@author: cosmo
"""

from qiskit import QuantumCircuit
import qiskit as qk

class Algorithm:
    #handles execution
    
    def __init__(self, shots, backend):
        self.shots = shots
        self.backend = backend
    
    def measure(self, qc, qb, cb):
        if isinstance(qb,list) and isinstance(cb,list):
            assert len(qb) == len(cb)
            #print("YES")
            for q,c in zip(qb,cb):
                
                qc.measure(q,c)
        else:
            #print(type(qc))
            qc.measure(qb,cb)
        
        job = qk.execute(qc, 
                        backend = self.backend, 
                        shots=self.shots)
        """
                        optimization_level=self.optimization_level,
                        noise_model=self.noise_model,
                        coupling_map=self.coupling_map,
                        basis_gates=self.basis_gates,
                        initial_layout=self.layout,
                        seed_transpiler=self.seed,
                        seed_simulator=self.seed)"""
    
        return job.result().get_counts(qc)
        
    def expectation(self, measurement):
        #get expectation value
        
        exp_value = 0
        for (state,value) in measurement.items():
            sign = 1 if state.count('1') % 2 == 0 else -1
            exp_value += sign * value/self.shots
        return exp_value
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        