
#include <mat.h>
#include "newton_cd_loop.h"
#include "newton_cd_loop_mat.h"

using Eigen::Map;
using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using std::complex;

MatrixXd GetMatrix(const mxArray* array, int m, int n) {
  CHECK(mxGetM(array) == m);
  CHECK(mxGetN(array) == n);
  return Map<MatrixXd>(mxGetPr(array), m, n);
}

MatrixXcd GetComplexMatrix(const mxArray* array, int m, int n) {
  CHECK(mxGetM(array) == m);
  CHECK(mxGetN(array) == n);
  double* real = mxGetPr(array);
  double* imag = mxGetPi(array);

  int p = 0;
  MatrixXcd B(m, n);
  for (int j = 0; j < n; j++) {
    for (int i = 0; i < m; i++, p++) {
      B(i,j) = complex<double>(real[p], imag[p]);
    }
  }

  return B;
}

const mxArray* _v;
#define GET_PARAM(pa, name, x, f) \
  _v = mxGetField(pa, 0, name);   \
  if (_v) x = f(_v)

void GetParams(const mxArray* a, NewtonCDParams* p) {
  GET_PARAM(a, "verbose", p->verbose, (int)mxGetScalar);
  GET_PARAM(a, "cd_max_iters", p->cd_max_iters, (int)mxGetScalar);
  GET_PARAM(a, "cd_tol", p->cd_tol, mxGetScalar);
  GET_PARAM(a, "cd_diagonal_only", p->cd_diagonal_only, (bool)mxGetScalar);
}
