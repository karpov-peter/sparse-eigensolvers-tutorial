// Example is based on https://github.com/opencollab/arpack-ng/blob/master/TESTS/icb_arpack_cpp.cpp

/*
 * This example demonstrates the use of C++ bindings to call arpack.
 *
 * Use arpack as you would have normally done, but, use [ae]upd instead
 * of *[ae]upd_. The main advantage is that compiler checks the argument types
 * and the correct function is called based on the type (float vs double vs
 * complex). Note: to debug arpack, call debug_c. This is a test program to
 * solve for the 9 eigenvalues of A*x = lambda*x where A is the diagonal
 * matrix with entries 1000, 999, ... , 2, 1 on the diagonal.
 */

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "arpack.hpp"

void diagonal_matrix_vector_product(const double* x, double* y) {
  for (int i = 0; i < 1000; ++i) {
    y[i] = static_cast<double>(i + 1) * x[i];
  }
}

void diagonal_matrix_vector_product(const std::complex<double>* x, std::complex<double>* y) {
  for (int i = 0; i < 1000; ++i) {
    // Use complex matrix (i, -i) instead of (i, i): this way "largest_magnitude"
    // and "largest_imaginary" options produce different results that can be checked.
    y[i] = x[i] * std::complex<double>{double(i + 1), -double(i + 1)};
  }
}

int main() {

  // arpack without debug
  //real_symmetric_runner<double>(1.e-05, arpack::which::largest_algebraic);

  const int N      = 1000;
  const int nev    = 9;
  const int ncv    = 2 * nev + 1;
  const int ldv    = N;
  const int ldz    = N;
  const int lworkl = ncv * (ncv + 8);
  const int rvec   = 1;      // need eigenvectors

  const double tol   = 1e-6; // small tol => more stable checks after EV computation.
  const double sigma = 0.0;  // not referenced in this mode

  std::vector<double> resid(N);
  std::vector<double> V(ldv * ncv);
  std::vector<double> z(ldz * nev);
  std::vector<double> d(nev);
  std::vector<double> workd(3 * N);
  std::vector<double> workl(lworkl);
  std::vector<int> select(ncv); // since HOWMNY = 'A', only used as workspace here

  int iparam[11], ipntr[11];
  iparam[0] = 1;      // ishift
  iparam[2] = 10 * N; // on input: maxit; on output: actual iteration
  iparam[3] = 1;      // NB, only 1 allowed
  iparam[6] = 1;      // mode

  arpack::which const ritz_option = arpack::which::largest_algebraic;

  int info = 0, ido = 0;
  do {
    arpack::saupd(ido, arpack::bmat::identity, N,
                  ritz_option, nev, tol, resid.data(), ncv,
                  V.data(), ldv, iparam, ipntr, workd.data(),
                  workl.data(), lworkl, info);

    diagonal_matrix_vector_product(&(workd[ipntr[0] - 1]), &(workd[ipntr[1] - 1]));
  } while (ido == 1 || ido == -1);

  // check info and number of ev found by arpack.
  if (info < 0 || iparam[4] < nev) { /*arpack may succeed to compute more EV than expected*/
    std::cout << "ERROR in saupd: iparam[4] " << iparam[4] << ", nev " << nev
              << ", info " << info << std::endl;
    throw std::domain_error("Error inside ARPACK routines");
  }

  arpack::seupd(rvec, arpack::howmny::ritz_vectors, select.data(), d.data(),
                z.data(), ldz, sigma, arpack::bmat::identity, N,
                ritz_option, nev, tol, resid.data(), ncv,
                V.data(), ldv, iparam, ipntr, workd.data(),
                workl.data(), lworkl, info);
  if (info < 0) throw std::runtime_error("Error in seupd, info " + std::to_string(info));

  for (int i = 0; i < nev; ++i) {
    double val = d[i];
    double ref = (N-(nev-1)+i);
    double eps = std::fabs(val - ref);
    std::cout << val << " - " << ref << " = " << eps << std::endl;

    /*eigen value order: smallest -> biggest*/
    if (eps > 1e-05) throw std::domain_error("Correct eigenvalues not computed");
  }
  std::cout << "Done" << std::endl;

  return 0;
}
