# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 16:45:57 2022

@author: cosmo
"""

a = random.uniform(0,2*np.pi)
b = random.uniform(0,2*np.pi)
desired_vector = [
    0,
    np.cos(a),
    np.sin(a)*np.sin(b),
    0,
    np.sin(a)*np.cos(b),
    0,
    0,
    0]


theta_range = np.linspace(0, 2 * np.pi, 128)

circuit = wavefunction.bind_parameters({alpha: a, beta: b})


print(state_fidelity(desired_vector,Statevector(circuit)))
circuit.measure([0,1,2],[0,1,2])
circuit.draw()


backend = Aer.get_backend('aer_simulator_statevector')
qc_compiled = transpile(circuit, backend)
job_sim = backend.run(qc_compiled, shots=1024)
result_sim = job_sim.result()
counts = result_sim.get_counts(qc_compiled)
print(counts)