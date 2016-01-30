
#ifndef NEWTON_CD_LOOP_MEX_H
#define NEWTON_CD_LOOP_MEX_H

MatrixXd GetMatrix(const mxArray* array, int m, int n);
MatrixXcd GetComplexMatrix(const mxArray* array, int m, int n);
void GetParams(const mxArray* array, NewtonCDParams* params);

#endif  // NEWTON_CD_LOOP_MEX_H
