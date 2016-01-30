clear params;
params.max_iters = 1000;
params.sigma = 1e-4;
params.beta = 0.5;
params.verbose = 1;
params.max_eig_limit = 0;
params.cd_tol = 1e-2;
params.cd_theta_eps = 1e-2;
params.tol = 1e-10;
params.ls_max_iters = 20;

lqrsp = @lqrsp_newton;
run_wac('newton', lqrsp, params, offdiag);
