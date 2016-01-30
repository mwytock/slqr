
clear params;
params.max_iters = 5000;
params.verbose = 1;
params.max_eig_limit = 0;
params.beta = 0.5;
params.sigma = 1e-4;
params.tol = 1e-9;
params.max_ls_iters = 20;
params.acceleration = false;

lqrsp = @lqrsp_fista;
run_springs('ista', lqrsp, params, offdiag);
