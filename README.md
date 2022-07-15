# VQE

The eigensolver is implemented in solver.py. It receives:

- an ansatz, created in ansatz.py. When it's called through __call__ it returns a quantum circuit containing a parametrized wavefunction. 
  It can either implement a generic UCCSD or the quarkonium optimized ansatz (if quark is True). The quarkonium one is all implemented with 2 parameters in UCCSD(), while 
  the general ones is constructed for a generic number of orbitals/fermions (hence parameters) by calling pauli_sum_op_to_exp_op_circuit(), which converts a pauli string 
  (obtained from mapping a combination of construction/destruption operators) into a circuit representing e^(pauli_string)
  This ansatz shape will be called from the vqe with different parameters each time
 
 - an hamiltonian, made of Pauli objects (basicly PauliOp objects from qiskit). It receives an array of coefficients, a number of fermions and orbitals (1 and 3 in the
   quarkonium case) and returns a list of Pauli. It's in hamiltonian.py
   
 - an optimizer, from optimizer.py, with the specs for optimization (the kind of optimizer, number of iterations, tol...)
 
 - an algorithm, from algorithm.py, with the specifics for the measuremnts (n. shots, if to use a simulator or a real QC, if include a noise model...)
 
 The vqe objects has several methods, mainly:
 
 - vqe_expval, that builds the ansatz with certain parameters, attaches one piece of hamiltonian at a time, measures the circuit through the algorithm and computes
   the expectation value through the pauli_string method (in hamiltonian.py). Then returns the total energy of the hamiltonian for that given ansatz parameters.
   
 - optimize_parameters optimizes with a classical optimizer the loss_function vqe_expval, so it finds the optimal parameters to return the smaller energy.
 
 - expval_excited_states is similar to vqe_expval (so it can be passed as loss function to an optimizer) but considers excited states. After finding the ground state 
   optimizing the parameters through vqe_expval, one needs to save the state via save_state(), so that the excited states will be orthogonal. The reason for that specific 
   return value can be found in 	arXiv:1805.08138
   
 - fold adds copies of (U_dagger U), where U is the complete quantum circuit, in order to increase the noise level (it has an effect only with a noisy device). This 
   allows to compute the ground state for different noise levels and extrapolate the noiseless limit. 
   Reference: T. Giurgica-Tiron, Y. Hindy, R. LaRose, A. Mari and W. J. Zeng, "Digital zero noise extrapolation for quantum error mitigation," 
   2020 IEEE International Conference on Quantum Computing and Engineering (QCE), 2020, pp. 306-316, doi: 10.1109/QCE49297.2020.00045.
   

quarkonium.py contains the spcifics for the quarkonium model hamiltonian, while mapping.py converts construction and destruction operators into pauli operators
through Jordan-Wigner mapping.

'main_excited' and 'main_noise' are the main files, respectively for finding ground + excited states and for finding the ground state with noise extrapoation.
In the mains, algorithm, optimizer, ansatz and hamiltonian objects need to be instantiated and passed to the eigensolver object. You can choose to use a real QC
or a simulator in the algorithm (IBMQ takes a lot more time, several hours at least to find the ground state). The ansatz can be a generic UCC (around 30 of circuit
depth for 3 orbitals and 1 fermion) or the quarkonium specific one that takes less time (set quark=True in the ansatz).  
