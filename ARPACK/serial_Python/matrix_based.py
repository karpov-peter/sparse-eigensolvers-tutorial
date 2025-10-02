import numpy as np
from scipy.sparse.linalg import eigsh
from scipy.sparse import diags

# Problem parameters
N = 1000        # matrix size
nev = 9         # number of requested eigenpairs
ncv = 2*nev + 1 # Krylov subspace size. ncv should be greater than nev. The recommended choice is ncv ~ 2*nev
tol = 1e-6      # relative accuracy for eigenvalues
maxiter = 10*N  # maximum number of iterations

# --- Sparse matrix (recommended) ---
# A = diag(1, 2, ..., N)
A = diags(np.arange(1, N+1, dtype=np.float64), offsets=0, format='csr')

# Largest algebraic value eigenpairs
vals, vecs = eigsh(A, k=nev, which='LA', ncv=ncv, tol=tol, maxiter=maxiter)

print("calculated eigenvalue - analytic eigenvalue = difference")

# Check correctness
ok = True
for i in range(nev):
    ref = float(N-(nev-1)+i)
    eps = abs(vals[i] - ref)
    print(f"{vals[i]:.6f} - {ref:.6f} = {eps:.6f}")
    if eps > 1e-5:
        print(f"Eigenvalue {i+1} does not match: {vals[i]} vs {ref}")
        ok = False

print("Done" if ok else "Failed")