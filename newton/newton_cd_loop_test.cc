
#include <assert.h>
#include <mat.h>
#include <sys/time.h>
#include <iostream>
#include <Eigen/Dense>
#include "newton_cd_loop.h"
#include "newton_cd_loop_mat.h"

using Eigen::Map;
using Eigen::MatrixXcd;
using Eigen::MatrixXd;

MATFile* file;
MatrixXd B, R, K, L, P, D, Lambda, active;
MatrixXcd U, V, X;
NewtonCDParams params;

double WallTime() {
  struct timeval time;
  CHECK(!gettimeofday(&time,NULL));
  return (double)time.tv_sec + (double)time.tv_usec * 1e-6;
}

mxArray* GetVariable(const char* name) {
  mxArray* array = matGetVariable(file, name);
  CHECK(array);
  return array;
}

void LoadInput() {
  // Set up sizes
  int m = mxGetM(GetVariable("K"));
  int n = mxGetN(GetVariable("K"));
  int r = mxGetN(GetVariable("X"));
  int k = mxGetM(GetVariable("active"));

  // Load variables
  B = GetMatrix(GetVariable("B"), n, m);
  R = GetMatrix(GetVariable("R"), m, m);
  K = GetMatrix(GetVariable("K"), m, n);
  L = GetMatrix(GetVariable("L"), n, n);
  P = GetMatrix(GetVariable("P"), n ,n);
  U = GetComplexMatrix(GetVariable("U"), n, n);
  V = GetComplexMatrix(GetVariable("V"), n, n);
  X = GetComplexMatrix(GetVariable("X"), n, r);
  Lambda = GetMatrix(GetVariable("Lambda"), m, n);
  active = GetMatrix(GetVariable("active"), k, 2);
  GetParams(GetVariable("params"), &params);
}

void CheckResults() {
  // TODO(mwytock): Implement this
}

void RunTest() {
  D = MatrixXd::Zero(K.rows(), K.cols());
  printf("Running NewtonCD (verbose = %d)\n", params.verbose);
  double start = WallTime();
  NewtonCD(B,R,K,L,P,U,V,X,Lambda,active,params,D);
  printf("Completed in %f seconds\n", WallTime()-start);
}

int main(int argc, char* argv[]) {
  if (argc != 2) {
    fprintf(stderr, "Usage: ./newton_cd_loop_test <input.mat>\n");
    return 1;
  }

  file = matOpen(argv[1], "r");
  CHECK(file);
  LoadInput();
  RunTest();
  if (matGetVariable(file, "D")) {
    CheckResults();
  }
  printf("PASSED\n");
}
