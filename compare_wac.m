% Compare algorithms on a simple mass spring example

load lqrsp/examples/neData;
W = B1*B1';
B = B2;

K0 = -lqr(A,B,Q,R);
lambda = 0.1;

% ADMM
clear params;
params.max_iters = 300;
params.rho = 100;
[K1,h1] = lqrsp_admm(A,B1,B,Q,R,-K0,lambda,params);

% FISTA
clear params;
params.max_iters = 300;
params.step = 1e-2;
params.verbose = 1;
[K2,h2] = lqrsp_fista(A,B,W,Q,R,lambda,K0,params);

% Newton method w/ coordinate descent
clear params;
params.max_iters = 20;
params.argmin_dir = @newton_cd;
params.sigma = 1e-4;
params.beta = 0.5;
params.verbose = 1;
params.proj_orthant = 0;
params.cd_max_iters = 20;
params.cd_tol = 1e-6;
params.cd_theta_eps = 1e-2;
[K4,h4] = lqrsp_newton(A,B,W,Q,R,lambda*ones(size(K0)),K0,params);

min_x = min([h1.objval h2.objval h4.objval]);
figure;
semilogy(h1.objval - min_x, 'LineWidth', 1);
hold on;
semilogy(h2.objval - min_x, 'Color', [0 0.5 0], 'LineWidth', 1);
semilogy(h4.objval - min_x, 'r', 'LineWidth', 1);
legend('ADMM', 'FISTA', 'Newton-CD');
prepare_figure('objval_wac.pdf', [8 6], 'Iteration', 'f - f^*');

figure;
semilogy(h1.time, h1.objval - min_x, 'LineWidth', 1);
hold on;
semilogy(h2.time, h2.objval - min_x, 'Color', [0 0.5 0], 'LineWidth', 1);
semilogy(h4.time, h4.objval - min_x, 'r', 'LineWidth', 1);
legend('ADMM', 'FISTA', 'Newton-CD');
prepare_figure('objval_wac_time.pdf', [8 6], 'Time (seconds)', 'f - f^*');
