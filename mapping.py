# -*- coding: utf-8 -*-
"""
Created on Tue Jun 21 16:21:50 2022

@author: cosmo
"""

from qiskit_nature.operators.second_quantization import FermionicOp
from qiskit_nature.mappers.second_quantization import JordanWignerMapper


def jw_mapping(n_qubits, qubit, coeff=1, dagger=True): #return a or a_dag
    label = ""
    for item in range(n_qubits):        
        if item == qubit:
            if dagger==True:
                label += "+"
            else:
                label += "-"
        else:
            label += "I"
    return FermionicOp((label, coeff), display_format="dense")

def one_body(n_qubits, a, i, coeff=1): #a^a_dag @ a_i
    op = jw_mapping(n_qubits, a, coeff, dagger=True) @ jw_mapping(n_qubits, i, dagger=False)
    op = JordanWignerMapper().map(op)
    return op
    
    
def two_body(n_qubits, a, b, i, j, coeff=1): #a_a+ @ a_b+ @ a_i @ a_j
    op = jw_mapping(n_qubits, a, coeff, dagger=True) @ jw_mapping(n_qubits, b, dagger=True) @ jw_mapping(n_qubits, i, dagger=False) @ jw_mapping(n_qubits, j, dagger=False)
    op = JordanWignerMapper().map(op)
    return op

                
                