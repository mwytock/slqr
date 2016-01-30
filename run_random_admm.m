% Script for running ADMM on random problem

clear params;
params.max_iters = 300;
params.rho = 100;
params.verbose = 1;
params.max_eig_limit = -1e-1;

lambda = 10;
n = 50;
[A,B,W,Q,R] = problem_random(n);
K0 = zeros(size(B'));
Lambda = lambda*ones(size(K));

[K1,h1] = lqrsp_admm(A,B,W,Q,R,-K0,Lambda,params);
