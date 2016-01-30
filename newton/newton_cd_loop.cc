// C++ implementation of the Newton-CD loop using Eigen
//
// Author: Matt Wytock (mwytock@cs.cmu.edu)

#include <complex>
#include <vector>
#include <iostream>
#include <Eigen/Dense>
#include "newton_cd_loop.h"

using Eigen::Map;
using Eigen::MatrixXcd;
using Eigen::MatrixXd;
using std::complex;
using std::vector;

double SoftThreshold(double a, double kappa) {
  return fmax(0, a-kappa) - fmax(0, -a-kappa);
}

class IterationState {
public:
  IterationState(const MatrixXd& B,
                 const MatrixXd& R_,
                 const MatrixXd& L_,
                 const MatrixXd& E,
                 const MatrixXcd& U,
                 const MatrixXcd& V,
                 const MatrixXcd& X)
    : R(R_),
      L(L_),
      ETL(E.transpose()*L),
      rank(X.cols()) {

    int m = B.cols();
    int n = B.rows();
    RD = MatrixXd::Zero(m,n);

    rank = X.cols();

    MatrixXcd Xk(X.rows(), X.rows());
    XB.resize(rank);
    XL.resize(rank);
    ETXB.resize(rank);
    ETXL.resize(rank);

    XLDT.resize(rank);
    DXL.resize(rank);
    XBD.resize(rank);
    DXB.resize(rank);

    for (int k = 0; k < rank; k++) {
      Xk = U*X.col(k).asDiagonal()*V;

      XB[k] = Xk*B;
      XL[k] = Xk*L;
      ETXB[k] = E.transpose()*XB[k];
      ETXL[k] = E.transpose()*XL[k];

      XLDT[k] = MatrixXcd::Zero(n,m);
      DXL[k] = MatrixXcd::Zero(m,n);
      XBD[k] = MatrixXcd::Zero(n,n);
      DXB[k] = MatrixXcd::Zero(m,m);
    }
  }

  void Update(double i, double j, double mu) {
    RD.col(j) += mu*R.col(i);
    for (int k = 0; k < rank; k++) {
      XLDT[k].col(i) += mu*XL[k].col(j);
      DXL[k].row(i)  += mu*XL[k].row(j);
      XBD[k].col(j)  += mu*XB[k].col(i);
      DXB[k].row(i)  += mu*XB[k].row(j);
    }
  }

  double ComputeSecondDifferential(const MatrixXd& D) {
    complex<double> d = (L*D.transpose()*R*D).trace();
    for (int k = 0; k < rank; k++) {
      complex<double> dk1 = (XL[k].transpose()*D.transpose()*ETXB[k]*D).trace();
      complex<double> dk2 = (ETXL[k].transpose()*D*XB[k]*D).trace();
      d -= 2. *(dk1 + dk2);
    }
    CHECK(fabs(d.imag()/d.real()) < 1e-8);
    return d.real();
  }

  // Compute the quadratic term, a
  double quadratic(int i, int j) const {
    complex<double> a = 2*R(i,i)*L(j,j);
    for (int k = 0; k < rank; k++) {
      a -= 4.*(ETXB[k](i,i)*XL[k](j,j) + ETXL[k](i,j)*XB[k](j,i));
    }
    CHECK(fabs(a.imag()/a.real()) < 1e-8);
    return a.real();
  }

  // Compute the linear term, b
  double linear(int i, int j, bool diagonal) const {
    complex<double> b = 2*ETL(i,j);
    if (!diagonal) {
      b += 2*RD.row(i)*L.col(j);
      for (int k = 0; k < rank; k++) {
        // NOTE(mwytock): Templating weirdness here!
        complex<double> bk1 = ETXB[k].row(i)*XLDT[k].transpose().col(j);
        complex<double> bk2 = ETXB[k].transpose().row(i)*DXL[k].col(j);
        complex<double> bk3 =  XBD[k].row(j)*ETXL[k].transpose().col(i);
        complex<double> bk4 = ETXL[k].transpose().row(j)*DXB[k].col(i);
        b -= 2.*(bk1 + bk2 + bk3 + bk4);
      }
    }
    CHECK(fabs(b.imag()/b.real()) < 1e-8);
    return b.real();
  }


private:
  MatrixXd R;
  MatrixXd L;
  MatrixXd ETL;

  // Rank of the Theta decomposition
  int rank;

  // Precomputation of various things involving the low rank Theta
  vector<MatrixXcd> XB;
  vector<MatrixXcd> XL;
  vector<MatrixXcd> ETXB;
  vector<MatrixXcd> ETXL;

  // Products involving D
  MatrixXd RD;
  vector<MatrixXcd> XLDT;
  vector<MatrixXcd> DXL;
  vector<MatrixXcd> XBD;
  vector<MatrixXcd> DXB;

  // Use diagonal Hessian approximation
  bool diagonal;
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
              MatrixXd& Ddiag) {
  bool indefinite = false;
  if (params.verbose >= 2)
    PRINTF("Eigen instruction sets: %s\n", Eigen::SimdInstructionSetsInUse());

  MatrixXd E = P*B + K.transpose()*R;
  IterationState state(B,R,L,E,U,V,X);

  for (int t = 0; t < params.cd_max_iters; t++) {
    MatrixXd Dold = D;
    for (int k = 0; k < active.rows(); k++) {
      int i = active(k,0);
      int j = active(k,1);

      double a = state.quadratic(i,j);
      // TODO(mwytock): Handle negative a case correctly
      if (a < 0)
        continue;
      
      if (!params.cd_diagonal_only) {
        double b = state.linear(i,j, false);
        double c = K(i,j) + D(i,j);
        double mu = -c + SoftThreshold(c - b/a, Lambda(i,j)/a);                
        D(i,j) = D(i,j) + mu;
        state.Update(i, j, mu);
      }

      if (t == 0) {
        double b = state.linear(i,j, true);
        double c = K(i,j);
        double mu = -c + SoftThreshold(c - b/a, Lambda(i,j)/a);
        Ddiag(i,j) = mu;
      }
    }

    double normD = D.norm();
    double diffD = (D - Dold).norm();
    if (params.verbose >= 2)
      PRINTF("CD %-3d %e %e\n", t, diffD, params.cd_tol*normD);
    if (diffD < params.cd_tol*normD)
      break;
  }
}
