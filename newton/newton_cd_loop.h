#ifndef NEWTON_CD_LOOP_H
#define NEWTON_CD_LOOP_H

#if (defined(LANG_M) || defined(MATLAB_MEX_FILE)) && 0
#include <mex.h>
#define PRINTF mexPrintf
#define CHECK(x) mxAssert(x, "Check failed")
#else
#include <stdio.h>
#define PRINTF printf
#define CHECK assert
#endif
#include <Eigen/Dense>

using Eigen::MatrixXd;
using Eigen::MatrixXcd;

struct NewtonCDParams {
NewtonCDParams()
: verbose(0),
    cd_max_iters(20),
    cd_tol(1e-6),
    cd_diagonal_only(false) {}

  // Verbosity level
  int verbose;

  // Stopping criteria
  int cd_max_iters;
  double cd_tol;

  // Do not compute normal direction, just diagonal
  bool cd_diagonal_only;
};

void NewtonCD(const MatrixXd& B,
              const MatrixXd& R,
              const MatrixXd& K,
              const MatrixXd& L,
              const MatrixXd& P,
              const MatrixXcd& U,
              const MatrixXcd& V,
              const MatrixXcd& X,
              const MatrixXd& Lambda,
              const MatrixXd& active,
              const NewtonCDParams& params,
              MatrixXd& D,
              MatrixXd& Ddiag);

#endif  // NEWTON_CD_LOOP_H
