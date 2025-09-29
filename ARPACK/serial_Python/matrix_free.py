import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

# Problem parameters
N    = 1000
k    = 9            # number of eigenpairs
ncv  = 2*k + 1      # Arnoldi subspace size
tol  = 1e-6
maxiter = 10*N

# Define A via a matvec: A = diag(1, 2, ..., N)
def dMatVec(x):
    # x is shape (N,)
    i = np.arange(1, N+1, dtype=np.float64)  # 1..N
    return i * x

A = LinearOperator(shape=(N, N), matvec=dMatVec, dtype=np.float64)

# Compute k smallest algebraic eigenpairs
# (since A is symmetric/Hermitian, use eigsh)
vals, vecs = eigsh(A, k=k, which='SA', ncv=ncv, tol=tol, maxiter=maxiter, return_eigenvectors=True)

# Check against exact eigenvalues 1..k
ok = True
for i in range(k):
    val = vals[i]
    ref = float(i+1)
    eps = abs(val - ref)
    print(f"{val:.6f} - {ref:.6f} = {eps:.6f}")
    if eps > 1e-5:
        print(f"Eigenvalue {i+1} does not match: {val} vs {ref}")
        ok = False

print("Done" if ok else "Failed") 
