
% Sweep over 100 values
lambdas = logspace(2, 0, 20);
[A,B,W,Q,R,K_mask] = problem_springs(n);

% Start from LQR with one FISTA step
K0 = -lqr(A,B,Q,R);
clear fparams;
fparams.max_iters = 1;
fparams.verbose = 1;
fparams.max_eig_limit = 0;
fparams.beta = 0.5;
fparams.sigma = 1e-4;
fparams.tol = 1e-2;
fparams.max_ls_iters = 100;
fparams.allow_unstable = true;
fparams.acceleration = false;
fparams.max_time = Inf;
Lambda = (1-K_mask)*lambdas(1);
K1 = lqrsp_fista(A,B,W,Q,R,Lambda,K0,fparams);

clear params;
params.max_iters = 100;
params.sigma = 1e-4;
params.beta = 0.5;
params.verbose = 1;
params.max_eig_limit = 0;
params.cd_tol = 1e-2;
params.cd_theta_eps = 1e-2;
params.tol = 1e-4;
params.ls_max_iters = 20;
params.max_time = Inf;

% Run the sweep with warm starts
Ks = {};
h = {};
Kpolish = {};
hpolish = {};
K = K1;
for i=1:length(lambdas)
  fprintf('i=%d\n', i);

  Lambda = ones(size(K_mask))*lambdas(i);
  [Ks{i},h{i}] = lqrsp_newton(A,B,W,Q,R,Lambda,K,params);
  K = Ks{i};
  Lambda = (K==0)*1e6;
  
  [Kpolish{i},hpolish{i}] = lqrsp_newton(A,B,W,Q,R,Lambda,K,params);
end

suffix = num2str(n);
save(['springs_sweep_' suffix '.mat'], 'lambdas', 'Ks', 'h', 'Kpolish', 'hpolish');
