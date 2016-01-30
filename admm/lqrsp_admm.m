function [K,h] = lqrsp_admm(A,B,W,Q,R,Lambda,K0,params)
% Solve LQRSP using ADMM with Anderson-Moore method for F minimization step.
[m,n] = size(K0);

% Primal variables
X = K0;
Z = K0;

% Dual variable
U = zeros(m,n);

% absolute and relative tolerances for the stopping criterion of ADMM
eps_abs = params.eps_abs;
eps_rel = params.eps_rel;
rho = params.rho;

tic;
for t=1:params.max_iters
  [At,U_,S_,V_,nu] = lqrsp_eig_deflate(A, B*X, params.max_eig_limit);

  params.am_max_iters = t;
  X = lqram(At,B,W,Q,R,Z-U,X,params); 
  Zold = Z;
  Z = st(X + U, Lambda/rho);
  U = U + (X - Z);

  r_norm = norm(X - Z, 'fro');
  s_norm = norm(-rho*(Z - Zold), 'fro');

  max_eig = max(real(eig(At + B*Z)));
  normK = sum(sum(abs(Lambda.*abs(X))));

  %  if max_eig < 0
  %   L = lyap(At + B*Z, W);
  %   J = trace(L*(Q + Z'*R*Z)) + normK;
  % else
  %   J = inf;
  % end
  L = lyap(At + B*X, W);
  J = trace(L*(Q + X'*R*X)) + normK;

  h.norm(t) = normK;
  h.nu(t) = nu;
  h.nnz(t) = nnz(Z);
  h.time(t) = toc;
  h.objval(t) = J;
    
  fprintf('ADMM %-4d %d %f %f %6.1e %6.1e\n', t, nnz(Z), nu, J, r_norm, s_norm);

  eps_pri  = sqrt(n*m)*eps_abs + eps_rel*max(norm(X,'fro'), norm(Z,'fro'));
  eps_dual = sqrt(n*m)*eps_abs + eps_rel*norm(rho*U,'fro');
  if  r_norm < eps_pri && s_norm < eps_dual
    break
  end

  if h.time(t) > params.max_time 
    if params.verbose >= 1
      fprintf('Time limit exceeded\n');
    end
    break
  end
end

K = Z;