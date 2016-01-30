% Compare algorithms on a simple mass spring example

N = 5;
I = eye(N,N);
Z = zeros(N,N);
T = toeplitz([2 -1 zeros(1,N-2)]);
A = [Z I; -T Z];
B1 = [Z; I];
W = B1*B1';
B = [Z; I];
Q = eye(2*N);
R = 10*I;

K0 = zeros(size(B'));
lambda = 1;
Lambda = lambda*ones(size(K0));

% % Newton method w/ conjugate gradient
% clear params;
% params.max_iters = 30;
% params.max_cg_iters = 20;
% params.argmin_dir = @newton_cg;
% params.sigma = 1e-4;
% params.beta = 0.5;
% params.verbose = 1;
% params.proj_orthant = 1;
% [K1,h1] = lqrsp_newton(A,B,W,Q,R,lambda,K0,params);

% % Newton method w/ coordinate descent (explicit Hessian)
% clear params;
% params.max_iters = 30;
% params.argmin_dir = @newton_cd_explicit;
% params.sigma = 1e-4;
% params.beta = 0.5;
% params.verbose = 1;
% params.proj_orthant = 0;
% params.cd_tol = 1e-6;
% params.cd_max_iters = 20;
% [K2,h2] = lqrsp_newton(A,B,W,Q,R,lambda,K0,params);

% Newton method w/ coordinate descent
clear params;
params.max_iters = 20;
params.argmin_dir = @newton_cd;
params.sigma = 1e-4;
params.beta = 0.5;
params.verbose = 1;
params.max_eig_limit = -1e-1;
params.tol = 1e-4;
params.cd_max_iters = 10;
params.cd_tol = 0.05;
params.cd_theta_eps = 1e-2;
[K3,h3] = lqrsp_newton(A,B,W,Q,R,Lambda,K0,params);

% min_x = min([h1.objval h2.objval h3.objval]);
% figure;
% semilogy(h1.objval - min_x, 'LineWidth', 1);
% hold on;
% semilogy(h2.objval - min_x, 'Color', [0 0.5 0], 'LineWidth', 1);
% semilogy(h3.objval - min_x, 'r', 'LineWidth', 1);
% legend('Newton-CG', 'Newton-CD (explicit)', 'Newton-CD');
% prepare_figure('objval_springs_small.pdf', [8 6], 'Iteration', 'f - f^*');
