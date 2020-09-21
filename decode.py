###################################################
# Decoding black hole interior from the radiation #
# 9/19/2020 Isaac H. Kim                          #
# MIT License                                     #
###################################################
import numpy as np

# Defining Pauli matrices
X = np.array([[0, 1], [1, 0]])
Y = np.array([[0, 1j], [-1j, 0]])
Z = np.array([[1, 0], [0, -1]])
I = np.array([[1, 0], [0, 1]])


def pauli_string(n_q, Xstring, Zstring):
    """
    Generates a matrix representation of a Pauli string

    Args:
        n_q(int): Number of qubits
        Xstring(int): Bitstring that stores whether X is
                      present or not on the i'th qubit.
        Zstring(int): Bitstring that stores whether Z is
                      present or not on the i'th qubit.
    """
    out = 1

    for i in range(n_q):
        if Xstring % 2:
            if Zstring % 2:
                out = np.kron(out, Y)
            else:
                out = np.kron(out, X)
        elif Zstring % 2:
            out = np.kron(out, Z)
        else:
            out = np.kron(out, I)
        Xstring //= 2
        Zstring //= 2
    return out


def apply_pauli_string(n_q, Xstring, Zstring, A):
    """
    Apply a Pauli string to a matrix A.

    Args:
        n_q(int): Number of qubits
        Xstring, Zstring(int): Bit string that represents the
                               Pauli-X and Pauli-Z matrices.
        A(np.ndarray): Matrix
    """
    A_o = A
    for i in range(n_q):
        if Xstring % 2:
            if Zstring % 2:
                P_i = Y
                P_i_embedded = np.kron(np.kron(np.eye(2**i), P_i),
                                       np.eye(2**(n_q-i-1)))
                A_o = P_i_embedded @ A_o
            else:
                P_i = X
                P_i_embedded = np.kron(np.kron(np.eye(2**i), P_i),
                                       np.eye(2**(n_q-i-1)))
                A_o = P_i_embedded @ A_o
        elif Zstring % 2:
            P_i = Z
            P_i_embedded = np.kron(np.kron(np.eye(2**i), P_i),
                                   np.eye(2**(n_q-i-1)))
            A_o = P_i_embedded @ A_o
    return A_o


class SparseDensity():
    """
    Sparse density matrix format. Represents the density
    matrix by its "square root." Specifically, let ρ
    be a (n x n) density matrix. We store it as

    ρ = A A^{†},

    where A is a (n x k) matrix, where k<<n.
    """
    def __init__(self, data=None):
        self.A = data

    @classmethod
    def rand(self, n_q, n_t):
        """
        Generate a random state over n_q qubits and trace out the last
        n_t qubits.
        """
        psi = np.random.normal(size=2**n_q) * 1j + np.random.normal(size=2**n_q)
        psi = psi / np.sqrt(sum(abs(psi)**2))

        sd = SparseDensity()
        sd.import_from_vec(psi, n_t)
        return sd

    @property
    def trace(self):
        return np.real(np.trace(self.A.conjugate().transpose() @ self.A))

    @property
    def dim(self):
        return self.A.shape[0]

    def overlap(self):
        """
        Compute the overlap with the |0> state on the first qubit.
        """
        zerovec = np.array([1, 0])
        vec = np.kron(zerovec, np.eye(self.dim//2))
        A_1 = vec @ self.A
        return np.real(np.trace(A_1.conjugate().transpose() @ A_1))

    def import_from_vec(self, psi, n_q):
        """
        Given a vector psi, trace out the last n_q qubits.
        """
        if len(psi.shape) > 1:
            raise TypeError("psi must be a vector.")
        if (len(psi) % (2**n_q)) != 0:
            raise TypeError("The dimension of psi not divisible by 2^(n_q).")

        dim = len(psi)
        self.A = np.reshape(psi, [dim//(2**n_q), 2**n_q])

    def coeffs(self, Xstring, Zstring):
        """
        The overlap of the first qubit with the |0> state, after applying
        ρ → exp(iPθ) ρ exp(-iPθ),
        must be of the following form:
        f(θ) = a_1 cos^2(θ) + a_2 sin^2(θ) + a_3cos(θ)sin(θ),
        where
        a_1 = Tr[<0|ρ|0>],
        a_2 = Tr[<0|PρP|0>],
        a_3 = -2Im(Tr[<0|Pρ|0>]),
        where |0> is the state of the first qubit. Using the cyclicity
        of the trace, we can recast these into the following (more
        efficiently computable) form.

        a_1 = Tr[A_1^† A_1],
        a_2 = Tr[A_2^† A_2],
        a_3 = -2Im(Tr[A_1^† A_2]),

        where A_1 = <0|A and A_2 = <0|PA.
        """
        n_q = len(bin(self.dim))-3
        P = pauli_string(n_q, Xstring, Zstring)
        zerovec = np.array([1, 0])
        vec = np.kron(zerovec, np.eye(self.dim//2))
        A_1 = vec @ self.A
        A_2 = (vec @ P) @ self.A

        a_1 = np.real(np.trace(A_1.conjugate().transpose() @ A_1))
        a_2 = np.real(np.trace(A_2.conjugate().transpose() @ A_2))
        a_3 = -2*np.imag(np.trace(A_1.conjugate().transpose() @ A_2))

        return a_1, a_2, a_3

    def optimize_theta(self, Xstring, Zstring):
        """
        Find θ that maximizes the overlap of the first qubit with the |0>
        state under the map ρ → exp(iPθ) ρ exp(-iPθ).

        The functional form of the overlap is
        (a_1+a_2)/2 + ((a_1-a_2)*cos(2θ) + a_3*sin(2θ))/2.

        The optimum is either
        θ = atan(a_3/(a_1-a_2))/2 or θ = atan(a_3/(a_1-a_2))/2 + pi/2,
        whichever yields the larger overlap.
        """
        a_1, a_2, a_3 = self.coeffs(Xstring, Zstring)

        theta1 = np.arctan(a_3/(a_1-a_2))/2
        theta2 = theta1 + np.pi/2

        f1 = ((a_1-a_2)*np.cos(2*theta1) + a_3*np.sin(2*theta1))/2
        f2 = ((a_1-a_2)*np.cos(2*theta2) + a_3*np.sin(2*theta2))/2

        f1 += (a_1+a_2)/2
        f2 += (a_1+a_2)/2

        if f1 > f2:
            return theta1, f1
        else:
            return theta2, f2

    def find_optimal_pauli(self):
        """
        Find the optimal Pauli and the optimal theta that maximizes
        the overlap of the first qubit with |0>.
        """
        i_max, j_max, f_max = 0, 0, 0
        theta_max = 0
        for i in range(self.dim):
            for j in range(self.dim):
                theta, f = self.optimize_theta(i, j)
                if f > f_max:
                    i_max, j_max, f_max = i, j, f
                    theta_max = theta
        return i_max, j_max, theta_max, f_max

    def update(self, Xstring, Zstring, theta):
        """
        Update A → exp(iPθ)A = cos(θ)A + isin(θ)PA.
        """
        n_q = len(bin(self.dim))-3
        P = pauli_string(n_q, Xstring, Zstring)
        self.A = np.cos(theta)*self.A + 1j*np.sin(theta)*(P@self.A)

    def optimal_update(self):
        """
        Find the optimal Pauli and theta and update.
        """
        f_old = self.overlap()
        i, j, theta, f = self.find_optimal_pauli()
        wt = bin(i | j).count('1')
        self.update(i, j, theta)
        print("Overlap changed {} → {}".format(f_old, f))
        return theta, f, wt

    def decode(self, eps):
        """
        Keep on decoding until the overlap is larger than 1-eps.
        """
        while self.overlap() < (1-eps):
            a, b = self.optimal_update()

    def decode_cost(self, eps):
        """
        Cost of the decoding up to an error of eps, quantified
        in terms of the
        1) Number of PPR(Pauli Product Rotations)
        2) Total amount of angles
        """
        itt = 0
        theta_sum = 0.0
        theta_sum_sq = 0.0
        weight_sum = 0
        while self.overlap() < (1-eps):
            theta, f, wt = self.optimal_update()
            itt += 1
            weight_sum += wt
            theta_sum += abs(theta)
            theta_sum_sq += abs(theta)**2
        theta_rms = np.sqrt(theta_sum_sq)
        return itt, theta_sum, theta_rms, weight_sum
