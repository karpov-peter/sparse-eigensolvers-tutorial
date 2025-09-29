import numpy as np
from scipy.sparse import diags
from scipy.sparse.linalg import eigsh

# Problem parameters
N = 1000
k = 9
ncv = 2*k + 1
tol = 1e-6
maxiter = 10*N

# --- Sparse matrix (recommended) ---
# A = diag(1, 2, ..., N)
A = diags(np.arange(1, N+1, dtype=np.float64), offsets=0, format='csr')

# Smallest algebraic eigenpairs
vals, vecs = eigsh(A, k=k, which='SA', ncv=ncv, tol=tol, maxiter=maxiter)

# Check against 1..k
ok = True
for i in range(k):
    ref = float(i+1)
    eps = abs(vals[i] - ref)
    print(f"{vals[i]:.6f} - {ref:.6f} = {eps:.6f}")
    if eps > 1e-5:
        print(f"Eigenvalue {i+1} does not match: {vals[i]} vs {ref}")
        ok = False

print("Done" if ok else "Failed")