% Verify the correctness of lqrsp_admm()

load ../lqrsp/examples/neData.mat;
lambda = 0.1;
F0 = lqr(A,B2,Q,R);

clear params;
params.max_iters = 300;
params.rho = 100;
[F1,h] = lqrsp_admm(A,B1,B2,Q,R,F0,lambda,params);

clear options;
options.method = 'l1';
options.gamval = lambda;
options.rho = 100;
options.maxiter = 100;
options.blksize = [1 1];
options.reweightedIter = 1;
solpath = lqrsp(A,B1,B2,Q,R,options);

assert(norm(F1-solpath.F,'fro') < 1e-4);
