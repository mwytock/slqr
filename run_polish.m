data = load('tac/results/springs_newton_500.mat');
[A,B,W,Q,R] = problem_springs(500);

K0 = data.K{1};
Lambda = (K0==0)*1e6;

clear params;
params.max_iters = 1000;
params.sigma = 1e-4;
params.beta = 0.5;
params.verbose = 2;
params.max_eig_limit = 0;
params.tol = 1e-10;
params.ls_max_iters = 100;
params.max_time = 3600;
params.max_cg_iters = nnz(K0);
params.use_lyap = false;
[K_cg,h_cg] = lqrsp_newton_cg(A,B,W,Q,R,Lambda,K0,params);

params.use_lyap = true;
[K_cg2,h_cg2] = lqrsp_newton_cg(A,B,W,Q,R,Lambda,K0,params);

clear params;
params.max_iters = 1000;
params.sigma = 1e-4;
params.beta = 0.5;
params.verbose = 1;
params.max_eig_limit = 0;
params.cd_tol = 1e-2;
params.cd_theta_eps = 1e-2;
params.tol = 1e-10;
params.ls_max_iters = 100;
params.max_time = Inf;

[K_cd,h_cd] = lqrsp_newton(A,B,W,Q,R,Lambda,K0,params);

save('springs_polish', 'K_cd', 'h_cd', 'K_cg', 'h_cg', 'K_cg2', 'h_cg2');