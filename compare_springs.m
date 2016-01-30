% Compare algorithms on a simple mass spring example

N = 50;
I = eye(N,N);
Z = zeros(N,N);
T = toeplitz([2 -1 zeros(1,N-2)]);
A = [Z I; -T Z];
B1 = [Z; I];
W = B1*B1';
B = [Z; I];
Q = eye(2*N);
R = 10*I;

lambda = 1;

K0 = zeros(size(B'));
Lambda = lambda*ones(size(K0));

% ADMM
clear params;
params.max_iters = 300;
params.rho = 100;
[K1,h1] = lqrsp_admm(A,B,W,Q,R,-K0,lambda,params);

% Newton method w/ coordinate descent

% [A,B,W,Q,R] = problem_springs();

% clear params;
% params.max_iters = 30;
% params.argmin_dir = @newton_cd;
% params.sigma = 1e-4;
% params.beta = 0.5;
% params.verbose = 2;
% params.proj_orthant = 0;
% params.max_eig_limit = -1e-1;
% params.cd_max_iters = 10;
% params.cd_tol = 0.05;
% params.cd_theta_eps = 1e-2;
% params.tol = 1e-4;

% K0 = zeros(size(B'));
% lambdas = logspace(3,-2,100);
% for i=1:length(lambdas)
%   Lambda = lambdas(i)*ones(size(K0));
%   [K0,h] = lqrsp_newton(A,B,W,Q,R,Lambda,K0,params);
% end

% min_x = min([h1.objval h3.objval]);
% figure;
% semilogy(h1.objval - min_x, 'LineWidth', 1);
% hold on;
% semilogy(h3.objval - min_x, 'r', 'LineWidth', 1);
% legend('ADMM', 'Newton-CD');
% prepare_figure('objval_springs.pdf', [8 6], 'Iteration', 'f - f^*');

% figure;
% semilogy(h1.time, h1.objval - min_x, 'LineWidth', 1);
% hold on;
% semilogy(h3.time, h3.objval - min_x, 'r', 'LineWidth', 1);
% legend('ADMM', 'Newton-CD');
% prepare_figure('objval_springs_time.pdf', [8 6], 'Time (seconds)', 'f - f^*');
