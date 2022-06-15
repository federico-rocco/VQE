# -*- coding: utf-8 -*-
"""
Created on Tue May 31 12:45:35 2022

@author: cosmo
"""

import qiskit as qk
qr = qk.QuantumRegister(2)
qc = qk.QuantumCircuit(qr)


from qiskit.opflow import CX, I, H, Zero, Plus, X, Y, Z, StateFn
bell = CX @ (I ^ H) @ Zero

bell = Plus^H@Zero


circ = (Y^Z)@StateFn(qc)
print(circ)
