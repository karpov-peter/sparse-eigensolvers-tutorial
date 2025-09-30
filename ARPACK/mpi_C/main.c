// Example is based on https://github.com/opencollab/arpack-ng/blob/master/PARPACK/TESTS/MPI/icb_parpack_c.c

/*
 * This example demonstrates the use of ISO_C_BINDING to call arpack
 * (portability). IMPORTANT: MPI communicators MUST be passed from C to Fortran
 * using MPI_Comm_c2f.
 *
 * Just use arpack as you would have normally done, but, use *[ae]upd_c instead
 * of *[ae]upd_. The main advantage is that compiler checks (arguments) are
 * performed at build time. Note: to debug parpack, call debug_c.
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "parpack.h"

/* test program to solve for the 9 largest eigenvalues of
 * A*x = lambda*x where A is the diagonal matrix
 * with entries 1000, 999, ... , 2, 1 on the diagonal.
 * */

void dMatVec(int start, int N, const double* x, double* y) {
  int i;
  for (i = 0; i < N; ++i) y[i] = ((double)(start + i + 1)) * x[i];
};


int main() {
  MPI_Init(NULL, NULL);

    const int N      = 1000;
  const int nev    = 9;
  const int ncv    = 2 * nev + 1;
  const int ldv    = N;
  const int ldz    = N;
  const int lworkl = ncv * (ncv + 8);
  const int rvec   = 1;        // need eigenvectors

  const double tol   = 1e-6; // small tol => more stable checks after EV computation.
  const double sigma = 0;    // not referenced in this mode

  double resid[N];
  double V[ldv * ncv];
  double z[ldz * nev];
  double d[nev];
  double workd[3 * N];
  double workl[lworkl];
  int select[ncv];

  int iparam[11], ipntr[11];
  iparam[0] = 1;       // ishift
  iparam[2] = 10 * N;  // on input: maxit; on output: actual iteration
  iparam[3] = 1;       // NB, only 1 allowed
  iparam[6] = 1;       // mode

  char bmat[]   = "I";
  char which[]  = "SM";
  char howmny[] = "A";
  
  MPI_Fint MPI_COMM_WORLD_fortran = MPI_Comm_c2f(MPI_COMM_WORLD);

  /// Split problem across each process/////////////////////
  int rank, nprocs;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

  int N_local = N / nprocs;
  if (rank < N % nprocs) // spread the remaining on each process
    N_local = N_local + 1;
  
  int start = rank*(N / nprocs) + ((rank < N%nprocs) ? rank : N%nprocs);
  printf("rank: %d, start: %d, size: %d\n", rank, start, N_local);
  
  //////////////////////////////////////////////////////////
  int info = 0, ido = 0;
  do {
    pdsaupd_c(MPI_COMM_WORLD_fortran, &ido, bmat, N_local, which, nev, tol, resid, ncv, V, ldv, iparam, ipntr, 
              workd, workl, lworkl, &info);

    dMatVec(start, N_local, &(workd[ipntr[0] - 1]), &(workd[ipntr[1] - 1]));
  } while (ido == 1 || ido == -1);

  // check info and number of ev found by arpack.
  if (info < 0 || iparam[4] < nev) {
    printf("Error in saupd: iparam[4] %d, nev %d, info %d\n", iparam[4], nev, info);
    return 1;
  }

  pdseupd_c(MPI_COMM_WORLD_fortran, rvec, howmny, select, d, z, ldz, sigma, bmat, N_local, which, nev, tol,
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
    printf("rank %d : %f - %f = %f\n", rank, val, ref, eps);

    if (eps > 1.e-05) return 1;
  }
  
  fflush(stdout);
  MPI_Barrier(MPI_COMM_WORLD);
  
  MPI_Finalize();
  return 0;
}