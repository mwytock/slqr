
clear params;
params.max_iters = 30;
params.argmin_dir = @newton_cd;
params.sigma = 1e-4;
params.beta = 0.5;
params.verbose = 2;
params.max_eig_limit = -1e-1;
params.cd_tol = 0.05;
params.cd_theta_eps = 1e-2;
params.tol = 1e-4;

[A,B,W,Q,R] = problem_random(2000);
K0 = zeros(size(B'));
lambda = 4.75;
Lambda = lambda*ones(size(K0));
[K,h] = lqrsp_newton(A,B,W,Q,R,Lambda,K0,params);

save('random2000_single.mat', 'lambda', 'h');
