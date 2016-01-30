
clear params;
params.max_iters = 30;
params.argmin_dir = @newton_cd;
params.sigma = 1e-4;
params.beta = 0.5;
params.verbose = 1;
params.max_eig_limit = -1e-1;
params.cd_max_iters = 10;
params.cd_tol = 0.05;
params.cd_theta_eps = 1e-2;
params.tol = 1e-4;

problem = @() problem_random(2000);
lambdas = logspace(1,-3,100);
sweep = lqrsp_sweep(problem, @lqrsp_newton, params, lambdas);

save('random.mat', 'lambdas', 'sweep');
