import sys
import numpy as np
import qutip as qt
from math import pi
from itertools import *
import matplotlib.pyplot as plt

from qiskit import QuantumCircuit, execute
from qiskit import Aer, IBMQ, transpile
from qiskit.providers.ibmq.managed import IBMQJobManager

#if not IBMQ.active_account():
#    IBMQ.enable_account('YOUR_TOKEN')

##############################################################################

if __name__ == "__main__":

    # For each Pauli string of the form "XIXXY", of length n,
    # construct the tensor product of the corresponding
    # Pauli operators. These form a basis for 2^n x 2^n Hermitian matrices.
    def construct_pauli_basis(n):
        IXYZ = {"I": qt.identity(2),\
                "X": qt.sigmax(),\
                "Y": qt.sigmay(),\
                "Z": qt.sigmaz()}
        return dict([("".join(p),\
                      qt.tensor(*[IXYZ[s] for s in p]))\
                    for p in product(IXYZ.keys(), repeat=n)])

    # Find the components of an operator in the Pauli basis.
    def op_pauli(O, basis):
        return dict([(s, (o.dag()*O).tr()/np.sqrt(len(basis)))\
                        for s, o in basis.items()])

    # Recovers the original matrix from the Pauli components.
    def pauli_op(P, basis):
        return sum([P[s]*o for s, o in basis.items()])

    # Pretty print for the Pauli decomposition.
    def print_pauli(P):
        for s, v in P.items():
            if not np.isclose(v, 0):
                print("%s: %.2f" % (s, v))

    ##############################################################################

    samples = 60
    n_shots = 4096
    local = True
    backend_name = "ibmq_5_yorktown" # "ibmq_16_melbourne" 

    def help():
        print("Usage: python vqe.py --samples 60 --n_shots 4096 --quantum/--classical/--backend <name>")
        print("Make sure to fill in your token!")
        sys.exit()         

    for i, arg in enumerate(sys.argv):
        try:
            if arg == "--samples":
                samples = int(sys.argv[i+1])
            if arg == "--n_shots":
                n_shots = int(sys.argv[i+1])
            if arg == "--backend":
                local = False
                backend_name = arg
            if arg == "--quantum":
                local = False
                backend_name = "ibmq_vigo"
            if arg == "--classical":
                local = False
                backed_name = "ibmq_qasm_simulator"
            if arg == "--help":   
                raise   
        except:
            help()

    ##############################################################################

    print("our matrix:")
    H = qt.Qobj(np.array([[1,0,0,0],\
                          [0,0,-1,0],\
                          [0,-1,0,0],\
                          [0,0,0,1]]))
    H.dims = [[2]*2, [2]*2]
    print(H)
    Hl, Hv = H.eigenstates()

    # Let's make sure we know the right answer.
    print("\ncompare:")
    print("lowest eigval =  %.4f" % (Hl[0]))
    print("eigenvector = \n%s" % (Hv[0]))
    print()

    basis = construct_pauli_basis(2)
    H_pauli = op_pauli(H, basis)

    ##############################################################################

    def ansatz(theta):
        circ = QuantumCircuit(2)
        circ.h(0)
        circ.cx(0,1)
        circ.rx(theta, 0)
        return circ

    ##############################################################################

    if local:
        backend = Aer.get_backend('qasm_simulator')
    else:    
        provider = IBMQ.load_account()
        job_manager = IBMQJobManager()
        backend = provider.get_backend(backend_name)

    print("laying out circuits...")
    # Creates circuits: for each theta, we run through all the operators 
    # in the Pauli decomposition, and create a circuit for each one
    # that starts with the ansatz, then consists of pre-measurement
    # rotations on each of the qubits, and then a full measurement.
    circs = []
    circ_indices = []
    thetas = np.linspace(0, 2*np.pi, samples)
    for t in thetas:
        for n, v in H_pauli.items():
            if not np.isclose(v, 0):
                circ = ansatz(t)
                for i, o in enumerate(n):
                    if o == "X":
                        circ.ry(-pi/2, i)
                    elif o == "Y":
                        circ.rx(pi/2, i)
                circ.measure_all()
                circs.append(circ)
                circ_indices.append((t, n))

    print("evaluating...\n")
    # We optimize the circuits, and send them off to be evaluated.
    circs = transpile(circs, backend=backend)
    jb = execute(circs, backend, shots=n_shots) if local \
            else job_manager.run(circs, backend=backend, name='vqe', shots=n_shots) 
    jb_results = jb.result() if local else jb.results()

    # We run through the results. For each Pauli operator P in the
    # decomposition, we get the counts which, when converted into probabilities,
    # weight the eigenvalues, which are the product of the {-1, 1}'s of each
    # qubit. We weight each of these ...Z... expectation values by
    # the component of P in the Pauli decomposition, and sum to get the 
    # full expectation value.
    # We do this for each theta, and find the theta with the least <H>.
    expectation_values = []
    for t in thetas:
        values = []
        for n, v in H_pauli.items():
            if not np.isclose(v, 0):
                if n.count("I") == len(n):
                    values.append(H_pauli[n])
                else:
                    counts = jb_results.get_counts(circ_indices.index((t, n)))
                    values.append(\
                        H_pauli[n]*\
                            sum([(c/n_shots)*\
                                    np.prod([1 if s_=='1' \
                                                else -1 for s_ in s])\
                            for s, c in counts.items()])\
                        )
        e = sum(values)
        print("theta = %.4f, <H> = %.8f" % (t, e))
        expectation_values.append(e)
    min_result = np.argmin(expectation_values)
    final_theta = thetas[min_result]
    print("\nfinal_theta: %.4f, <H> = %.4f" % (final_theta, expectation_values[min_result]))

    # Use a statevector simulator to look at our final state.
    state_backend = Aer.get_backend('statevector_simulator')
    state = execute(ansatz(final_theta), backend=state_backend).result().get_statevector()
    print("final vector = \n%s" % qt.Qobj(state))

    # Plot!
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    ax.plot(thetas, expectation_values)
    ax.set_xlabel('θ', fontsize=14)
    ax.set_ylabel('<H>', fontsize=14)
    plt.show()