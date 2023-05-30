import numpy as np
from scipy.sparse import linalg 

m = 3
it = m + 1
A = np.array([[1, 7, 5], [8,2,0], [0,1,3]])
xp = np.array([1,1,1])
b = np.matmul(A,xp)
x0 = np.array([0,0,0])

C = np.array([[1, 0, 4, 0], 
              [0, 2, 2, 0],
              [4, 2, 1, 0],
              [0, 0, 0, 2]])

print("Test vectorroducts")



def arnoldi_iteration(A, b, n: int):
    """Compute a basis of the (n + 1)-Krylov subspace of the matrix A.

    This is the space spanned by the vectors {b, Ab, ..., A^n b}.

    Parameters
    ----------
    A : array_like
        An m × m array.
    b : array_like
        Initial vector (length m).
    n : int
        Dimension of the Krylov subspace, must be >= 1.
    
    Returns
    -------
    Q : numpy.array
        An m x (n + 1) array, where the columns are an orthonormal basis of the Krylov subspace.
    h : numpy.array
        An (n + 1) x n array. A on basis Q. It is upper Hessenberg.
    """
    eps = 1e-12
    h = np.zeros((n + 1, n))
    Q = np.zeros((A.shape[0], n + 1))
    # Normalize the input vector
    Q[:, 0] = b / np.linalg.norm(b, 2)  # Use it as the first Krylov vector
    for k in range(1, n + 1):
        v = np.dot(A, Q[:, k - 1])  # Generate a new candidate vector
        for j in range(k):  # Subtract the projections on previous vectors
            h[j, k - 1] = np.dot(Q[:, j].conj(), v)
            v = v - h[j, k - 1] * Q[:, j]
        h[k, k - 1] = np.linalg.norm(v, 2)
        if h[k, k - 1] > eps:  # Add the produced vector to the list, unless
            Q[:, k] = v / h[k, k - 1]
        else:  # If that happens, stop iterating.
            return Q, h
    return Q, h


def gmres(A, b, x0, tol=1e-6, max_iter=100):
    n = len(b)
    m = min(max_iter, n)  # Maximum number of iterations

    Q = np.zeros((n, m + 1))
    H = np.zeros((m + 1, m))
    c = np.zeros(m + 1)
    s = np.zeros(m + 1)

    r0 = b - np.dot(A, x0)
    beta = np.linalg.norm(r0)
    Q[:, 0] = r0 / beta

    for j in range(m):
        v = np.dot(A, Q[:, j])
        for i in range(j + 1):
            H[i, j] = np.dot(Q[:, i], v)
            v = v - H[i, j] * Q[:, i]

        H[j + 1, j] = np.linalg.norm(v)
        Q[:, j + 1] = v / H[j + 1, j]

        # Apply Givens rotation
        for k in range(j):
            temp = c[k] * H[k, j] + s[k] * H[k + 1, j]
            H[k + 1, j] = -s[k] * H[k, j] + c[k] * H[k + 1, j]
            H[k, j] = temp

        # Compute Givens rotation
        rho = np.linalg.norm([H[j, j], H[j + 1, j]])
        c[j] = H[j, j] / rho
        s[j] = H[j + 1, j] / rho
        H[j, j] = rho

        # Update residual norm
        beta = abs(s[j]) * beta
        if beta < tol:
            y = np.linalg.solve(H[:j+1, :j], beta * np.eye(j+1))
            return x0 + np.dot(Q[:, :j], y)

    y = np.linalg.lstsq(H[:m, :m], beta * np.eye(m))[0]
    return x0 + np.dot(Q[:, :m], y)

Q, h = arnoldi_iteration(A, b, m)
print("----------------ref---------")
print("cond number A", np.linalg.cond(A))
print(Q)
print(h)
print(h[0,1])
print(np.dot(Q[0], Q[1]))
print(np.dot(Q[0], Q[2]))
print(np.dot(Q[1], Q[2]))
print("----------------------------")

x, exitc = linalg.gmres(A, b, maxiter = 1)
#x = gmres(A, b, x0, max_iter = 10)

print("----------------ref---------")
print("Solution:", x)
print(exitc)
print("----------------------------")

#R = np.array([[9.746, 0.534522], [0, 6.04743]])
R = np.array([[1, 2.3, 0.782, 212], [0, 1, 0, 7], [0, 0, 3 , 1], [0 , 0 , 0 ,4]])
print(R)
g = np.array([1.59934, -0.125656, 5.234, 2.13123])
xa = np.linalg.solve(R, g)
b_test = np.matmul(R, xa)
print(xa)
print(b_test)

B = np.array([[1, 0, 4], 
              [0, 2, 2], 
              [4, 2, 1]])
C = np.array([[1,0,0], [0, 2, 0], [4, 2, 1]])
print(np.linalg.inv(C))
print("----------------ref CG---------")
b = np.matmul(B, xp)
print(b)
xc = linalg.cg(B, b, maxiter = 2)
print(xc)

y = np.matmul(np.linalg.inv(C),A)
print(y)

print("------------forwardSub + backwardSub--------")
xp = np.array([1.3, 1, 1, 9.2])
D = np.array([[1,0,0,0], [5,6,0,0], [9,10.4,11,0], [13, 14,15,16]])
k = np.matmul(D, xp)
print(k)

F = np.array([[1, 2.4, 3.4, 4], [0, 6, 7, 8], [0,0,11,12], [0,0,0,16]])
b = np.matmul(F, xp)
print(b)

#E = np.array([[1, 7, 5], [0, 2, 0],[0,0,3]])

print(np.linalg.solve(D, k))
print(np.linalg.solve(F, b))