
#include <mex.h>
#include <Eigen/Dense>
#include "newton_cd_loop.h"
#include "newton_cd_loop_mat.h"

using Eigen::Map;
using Eigen::MatrixXcd;
using Eigen::MatrixXd;

enum InputArguments {
  INPUT_B,
  INPUT_R,
  INPUT_K,
  INPUT_L,
  INPUT_P,
  INPUT_U,
  INPUT_V,
  INPUT_X,
  INPUT_LAMBDA,
  INPUT_ACTIVE,
  INPUT_PARAMS,
  NUM_INPUT
};

enum OutputArguments {
  OUTPUT_D,
  OUTPUT_D_DIAG,
  NUM_OUTPUT
};

// Helper conversion function
mxArray* GetMxArray(const MatrixXd& A) {
  mxArray* B = mxCreateDoubleMatrix(A.rows(), A.cols(), mxREAL);
  memcpy(mxGetPr(B), A.data(), A.rows()*A.cols()*mxGetElementSize(B));
  return B;
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {
  if (nrhs != NUM_INPUT)
    mexErrMsgIdAndTxt("NewtonCD:arguments", "Wrong number of input arguments.");
  if (nlhs > NUM_OUTPUT)
    mexErrMsgIdAndTxt("NewtonCD:arguments", "Wrong number of output arguments.");

  int m = mxGetM(prhs[INPUT_K]);
  int n = mxGetN(prhs[INPUT_K]);
  int r = mxGetN(prhs[INPUT_X]);
  int k = mxGetM(prhs[INPUT_ACTIVE]);

  MatrixXd B = GetMatrix(prhs[INPUT_B], n, m);
  MatrixXd R = GetMatrix(prhs[INPUT_R], m, m);
  MatrixXd K = GetMatrix(prhs[INPUT_K], m, n);
  MatrixXd L = GetMatrix(prhs[INPUT_L], n, n);
  MatrixXd P = GetMatrix(prhs[INPUT_P], n ,n);
  MatrixXcd U = GetComplexMatrix(prhs[INPUT_U], n, n);
  MatrixXcd V = GetComplexMatrix(prhs[INPUT_V], n, n);
  MatrixXcd X = GetComplexMatrix(prhs[INPUT_X], n, r);
  MatrixXd Lambda = GetMatrix(prhs[INPUT_LAMBDA], m, n);
  MatrixXd active = GetMatrix(prhs[INPUT_ACTIVE], k, 2);

  NewtonCDParams params;
  GetParams(prhs[INPUT_PARAMS], &params);

  MatrixXd D = MatrixXd::Zero(m,n);
  MatrixXd Ddiag = MatrixXd::Zero(m,n);
  NewtonCD(B,R,K,L,P,U,V,X,Lambda,active,params,D,Ddiag);
  if (nlhs > OUTPUT_D)
    plhs[OUTPUT_D] = GetMxArray(D);
  if (nlhs > OUTPUT_D_DIAG)
    plhs[OUTPUT_D_DIAG] = GetMxArray(Ddiag);
}
