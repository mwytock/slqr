% Script for running FISTA on random data

clear params;
params.max_iters = 5000;
params.verbose = 1;
params.max_eig_limit = 0; 
params.beta = 0.5;
params.sigma = 1e-4;
params.tol = 1e-6;
params.max_ls_iters = 20;
params.acceleration = true;

lqrsp = @lqrsp_fista;
run_wac('fista', lqrsp, params);
