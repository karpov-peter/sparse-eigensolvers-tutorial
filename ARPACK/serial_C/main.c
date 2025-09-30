// Example is based on https://github.com/opencollab/arpack-ng/blob/master/TESTS/icb_arpack_c.c

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "arpack.h"

/* test program to solve for the 9 smallest eigenvalues of
 * A*x = lambda*x where A is the diagonal matrix
 * with entries 1000, 999, ... , 2, 1 on the diagonal.
 * */

void dMatVec(const double* x, double* y) {
  int i;
  for (i = 0; i < 1000; ++i) y[i] = ((double)(i + 1)) * x[i];
};


int main() {
  
  const int N      = 1000;
  const int nev    = 9;
  const int ncv    = 2 * nev + 1; // Arnoldi subspace size. ncv should be greater than nev. The recommended choice is ncv ~ 2*nev
  const int ldv    = N;
  const int ldz    = N;
  const int lworkl = ncv * (ncv + 8);
  const int rvec   = 1;      // need eigenvectors

  const double tol = 1e-6; // small tol => more stable checks after EV computation.
  const double sigma = 0;      // not referenced in this mode
  
  double resid[N];
  double V[ldv * ncv];
  double z[ldz * nev];
  double d[nev]; // eigenvalues
  double workd[3 * N];
  double workl[lworkl];
  int select[ncv]; // since HOWMNY = 'A', only used as workspace here

  int iparam[11], ipntr[11];
  iparam[0] = 1;       // ishift
  iparam[2] = 10 * N;  // on input: maxit; on output: actual iteration
  iparam[3] = 1;       // NB, only 1 allowed
  iparam[6] = 1;       // mode

  char bmat[]   = "I"; // I -- standard eigenproblem; G -- generalized eigenproblem
  char which[]  = "SA"; // SA/LA -- smallest/largest algebraic value; SM/LM -- smallest/largest magnitude; BE -- both ends (nev/2 from each end)
  char howmny[] = "A"; // A -- all nev eigenvectors to be computed; S -- selected eigenvectors to be computed

  int info = 0, ido = 0;
  do {
    // "double symmetric Arnoldi update"
    dsaupd_c(&ido, bmat, N, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr,
             workd, workl, lworkl, &info);

    dMatVec(&(workd[ipntr[0] - 1]), &(workd[ipntr[1] - 1]));
  } while (ido == 1 || ido == -1);

  // check info and number of ev found by arpack.
  if (info < 0 || iparam[4] < nev) {
    printf("Error in saupd: iparam[4] %d, nev %d, info %d\n", iparam[4], nev, info);
    return 1;
  }

  // If dsaupd has converged, we retrieve the eigenvalues and compute the eigenvectors with dseupd.
  // "double symmetric eigenvectors update"
  dseupd_c(rvec, howmny, select, d, z, ldz, sigma, bmat, N, which, nev, tol,
           resid, ncv, V, ldv, iparam, ipntr, workd, workl, lworkl, &info);
  if (info < 0) {
    printf("Error in seupd: info %d\n", info);
    return 1;
  }

  int i;
  for (i = 0; i < nev; ++i) {
    double val = d[i];
    double ref = i+1;
    double eps = fabs(val - ref);
    printf("%f - %f = %f\n", val, ref, eps);

    /*eigen value order: smallest -> biggest*/
    if (eps > 1.e-05) 
      {
      printf("Eigenvalue %d does not match: %f vs %f\n", i, val, ref);
      return 1;
      }
  }

  printf("Done\n");
  return 0;
} 
