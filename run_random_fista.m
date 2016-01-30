% Script for running FISTA on random data

clear params;
params.max_iters = 300;
params.verbose = 1;
params.max_eig_limit = -1e-1;
params.beta = 0.5;
params.sigma = 1e-4;
params.tol = 1e-4;

% % Small examples, n=[50 100 500], whole regularization paths
% ns=[50,100,500];
% for i=1:length(ns)
%   problem = @() problem_random(ns(i));
%   lambdas = logspace(3,-3,100);
%   sweep = lqrsp_sweep(problem, @lqrsp_fista, params, lambdas);
%   save(['random_fista_' num2str(ns(i)) '.mat'], 'lambdas', 'sweep');
% end

lambda = 10;
n = 50;
[A,B,W,Q,R] = problem_random(n);
K0 = zeros(size(B'));
Lambda = lambda*ones(size(K));

[K,h] = lqrsp_fista(A,B,W,Q,R,Lambda,K0,params);
