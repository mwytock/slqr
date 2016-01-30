
K0 = -lqr(A,B,Q,R);
P = lyap((A + B*K0)', Q + K0'*R*K0);
J_lqr = trace(P*W)

clear params;
params.max_iters = 1000;
params.sigma = 1e-4;
params.beta = 0.5;
params.verbose = 1;
params.max_eig_limit = -1e-1;
params.cd_tol = 1e-2;
params.cd_theta_eps = 1e-2;
params.tol = 1e-6;
params.ls_max_iters = 20;
params.max_time = 100;

Lambda = (1-K_mask)*1e6;
K1 = lqrsp_newton(A,B,W,Q,R,Lambda,K0,params);
P = lyap((A + B*K1)', Q + K1'*R*K1);
J_local = trace(P*W)

loss = (J_local - J_lqr)/J_lqr