% Script for running Newton on springs data

clear params;
params.rho = 100;
params.eps_abs = 1e-8;
params.eps_rel = 1e-8;
params.max_iters = 10000;
params.sigma = 0.3;
params.beta = 0.5;
params.verbose = 1;
params.max_eig_limit = 0;
params.am_max_iters = 100;
params.am_ls_max_iters = 100;
params.am_tol = 1e-4;

lqrsp = @lqrsp_admm;
run_wac('admm', lqrsp, params, offdiag);
