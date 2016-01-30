function [K,h] = lqrsp_newton_cg(A,B,W,Q,R,Lambda,K0,params)
% Solve LQRSP using Newton-Lasso
K = K0;
J = Inf;

% Count the number of steps with nu < 0
tnu = 1;

tic;
for t=1:params.max_iters
  [At,U,S,V,nu] = lqrsp_eig_deflate(A, B*K, params.max_eig_limit);

  %L = lyap_fast(U, S, V, W);
  %P = lyap_fast(V.', S, U.', Q + K'*R*K);
  L = lyap(At + B*K, W);
  P = lyap((At + B*K)', Q + K'*R*K);
  normK = sum(sum(Lambda.*abs(K)));
  Jp = J;
  J = trace(L*(Q + K'*R*K)) + normK;

  % Compute the gradient
  G = 2*(R*K + B'*P)*L;
  subgrad = sum(sum(abs(~~K.*(G + sign(K).*Lambda) + ...
                        ~K.*max(abs(G) - Lambda, 0))));

  h.nnz(t) = nnz(K);
  h.nu(t) = nu;
  h.objval(t) = J;
  h.norm(t) = normK;
  h.time(t) = toc;

  if params.verbose >= 1
    fprintf('Newton %-3d %d %f %e\n', t, nnz(K), J, nu);
  end

  if h.time(t) > params.max_time
    if params.verbose >= 1
      fprintf('Time limit exceeded\n');
    end
    break
  end

  if subgrad < params.tol || subgrad < params.tol*normK || abs((Jp-J)/J) < 2.22e-16
    if params.verbose >= 1
      fprintf('Converged after %d iterations\n', t);
    end
    break;
  end

  params.cd_max_iters = ceil(tnu/3);
  params.cd_diagonal_only = false;
  if nu > 0
    params.cd_diagonal_only = true;
  else
    tnu = tnu + 1;
  end

  D = newton_cg(A,U,S,V,B,W,Q,R,Lambda,K,L,P,G,params);
  K = lqrsp_newton_ls(At,B,W,Q,R,Lambda,J,G,K,D,params);
end
