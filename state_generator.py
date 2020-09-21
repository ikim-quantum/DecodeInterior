###################################################
# Decoding black hole interior from the radiation #
# 9/19/2020 Isaac H. Kim                          #
# MIT License                                     #
###################################################
import numpy as np

X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, 1j], [-1j, 0]])
Z = np.array([[1, 0], [0, -1]])
I = np.array([[1, 0], [0, 1]])


def randU(n):
    """
    Generates a random (n x n) unitary matrix.

    Args:
        n(int): Dimension of the vector space

    Returns:
        np.ndarray: 4 x 4 unitary matrix
    """
    X = np.random.randn(n, n) + 1j * np.random.randn(n, n)
    Q, R = np.linalg.qr(X)
    R = np.diag(np.diag(R) / abs(np.diag(R)))
    U = Q @ R
    return U


class Gate():
    """
    Gate class
    Attrs:
        q1, q2(int): Qubit index
        U(np.ndarray): 4 x 4 unitary matrix.
    """
    def __init__(self, q1, q2, U):
        self.q1, self.q2, self.U = q1, q2, U

    @classmethod
    def rand(self, q1, q2):
        """
        Generates a random gate.

        Returns:
            Gate: Random gate
        """
        return Gate(q1, q2, randU(4))

    def __repr__(self):
        return "Gate over ({}, {})".format(self.q1, self.q2)


class Circuit():
    """
    Circuit class
    Attrs:
        circ(list(Gate)): A sequence of gates.
    """
    def __init__(self, circ):
        self.circ = circ

    @property
    def n_q(self):
        """
        Returns:
            Number of qubits in the circuit.
        """
        qs_max = [max(g.q1, g.q2) for g in self.circ]
        return max(qs_max)+1

    @classmethod
    def rand_1d(self, n_q, d):
        """
        Generates a random depth-d quantum circuit over n_q qubits.

        Args:
            n_q(int): Number of qubits
            d(int): Depth

        Returns:
            Circuit: Depth-d random 1D circuit over n_q qubits
        """
        circ = []
        pairs1 = [(i, (i + 1) % n_q) for i in range(n_q) if i % 2 == 0]
        pairs2 = [(i, (i + 1) % n_q) for i in range(n_q) if i % 2 == 1]
        for i in range(d):
            if i % 2 == 0:
                pairs = pairs1
            else:
                pairs = pairs2
            for pair in pairs:
                q1, q2 = pair
                circ.append(Gate.rand(q1, q2))
        return Circuit(circ)

    def run(self):
        return run(self)


def apply_gate(g, psi):
    """
    Apply gate g in on a state psi.

    Args:
        g(Gate): Gate
        psi(np.ndarray): A state vector

    Returns:
        np.ndarray: g @ psi
    """
    dim = len(psi)
    n = 0
    while dim > 1:
        if dim % 2 != 0:
            raise ValueError("Dimension of psi must be a power of two.")
        dim //= 2
        n += 1
    if g.q1 > n or g.q2 > n:
        raise ValueError("Qubit index out of range.")
    psi_tensor = np.reshape(psi, [2]*n)
    if g.q2 > 0:
        perm = [i for i in range(n)]
        perm[g.q1], perm[g.q2-1] = perm[g.q2-1], perm[g.q1]
        psi_tensor = np.transpose(psi_tensor, axes=perm)
        g_embedded = np.kron(np.kron(np.eye(2**(g.q2-1)), g.U),
                             np.eye(2**(n-g.q2-1)))
        psi_temp = g_embedded @ np.reshape(psi_tensor, 2**n)
        psi_tensor2 = np.reshape(psi_temp, [2]*n)
        psi_tensor2 = np.transpose(psi_tensor2, axes=perm)
        return np.reshape(psi_tensor2, 2**n)
    elif g.q1 > 0:
        perm1 = [i for i in range(n)]
        perm2 = [i for i in range(n)]
        perm1[g.q2], perm1[g.q1-1] = perm1[g.q1-1], perm1[g.q2]
        perm2[g.q1], perm2[g.q1-1] = perm2[g.q1-1], perm2[g.q1]
        psi_tensor = np.transpose(np.transpose(psi_tensor, axes=perm1),
                                  axes=perm2)
        g_embedded = np.kron(np.kron(np.eye(2**(g.q1-1)), g.U),
                             np.eye(2**(n-g.q1-1)))
        psi_temp = g_embedded @ np.reshape(psi_tensor, 2**n)
        psi_tensor2 = np.reshape(psi_temp, [2]*n)
        psi_tensor2 = np.transpose(np.transpose(psi_tensor2, axes=perm2),
                                   axes=perm1)
        return np.reshape(psi_tensor2, 2**n)


def run(c):
    """
    Run a Circuit instance, starting with the all-0 state.
    """
    psi = np.zeros(2**c.n_q)
    psi[0] = 1
    for g in c.circ:
        psi = apply_gate(g, psi)
    return psi
