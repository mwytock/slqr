% Script for running Newton on random data

clear params;
params.max_iters = 100;
params.argmin_dir = @newton_cd;
params.sigma = 1e-4;
params.beta = 0.5;
params.verbose = 1;
params.max_eig_limit = 0;;
params.tol = 1e-2;
params.ls_max_iters = 20;
params.cd_tol = 0.05;
params.cd_theta_eps = 1e-2;

% Small examples, whole regularization path, starting from middle
lambdas = logspace(-2,2,100);

ns=[500,100,50];
lqrsp = @lqrsp_newton;
for i=1:length(ns)
  problem = @() problem_random(ns(i));
  sweep = lqrsp_sweep(problem, lqrsp, params, lambdas);
  save(['random_newton_' num2str(ns(i)) '.mat'], 'lambdas', 'sweep');
end

% Sparse example w/ Polishing
% lambda = 10;
% n = 50;
% [A,B,W,Q,R] = problem_random(n);
% K0 = zeros(size(B'));
% Lambda = lambda*ones(size(K));
% [K,h] = lqrsp_newton(A,B,W,Q,R,Lambda,K0,params);

% Lambda = 1e6*(K == 0);
% K0 = K;
% [K,h] = lqrsp_newton(A,B,W,Q,R,Lambda,K0,params);

% % Specific examples from zero
% lambda = 1;
% n = 500;
% [A,B,W,Q,R] = problem_random(n);
% K0 = zeros(size(B'));
% Lambda = lambda*ones(size(K0));

% K = K0;
% while true
%   [At,U,S,V,nu] = lqrsp_eig_deflate(A, B*K, params.max_eig_limit);
%   fprintf('nu = %f\n', nu);
%   [K,h] = lqrsp_newton(At,B,W,Q,R,Lambda,K,params);

%   Lambda = 1e6*(K == 0);
%   [K,h] = lqrsp_newton(A,B,W,Q,R,Lambda,K,params);
% end
