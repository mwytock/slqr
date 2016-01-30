

load ../lqrsp/examples/neData.mat;
lambda = 0.1;
F0 = lqr(A,B2,Q,R);
clear params;
params.max_iters = 300;
params.step = 1e-2;
[F,h] = lqrsp_fista(A,B1,B2,Q,R,F0,lambda,params);
