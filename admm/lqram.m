function K = lqram(A,B,W,Q,R,U,K0,params)
[m,n] = size(K0);
% Solve L2 regularized LQR using Anderson-Moore method
%
% minimize f(K) = J(K) + (rho/2)||K - U||_F^2
%
rho = params.rho;

% check if F0 is a stabilizing feedback gain
max_eig  = max(real(eig(A + B*K0)));
if max_eig > 0
  error('Initial condition K0 is not stabilizing in Anderson-Moore method!')
end

IR = inv(R);
K = K0;

% controllability gramians and objective function phi
L = lyap(A + B*K, W);

f = trace(L*(Q + K'*R*K)) + (rho/2)*norm(K - U, 'fro')^2;
for k=1:params.am_max_iters
  % observability gramian
  P = lyap((A + B*K)', Q + K'*R*K);

  % one Sylvester equation for F
  Kbar = lyap(rho*IR, 2*L, IR*(2*B'*P*L - rho*U));

  % descent direction Kt
  Kt = Kbar - K;

  % gradient direction
  G = 2*(R*K + B'*P)*L + rho*(K - U);

  % check if Kt is a descent direction;
  if trace(Kt'*G) > 1e-10
    error('Kt is not a descent direction!');
  end

  normG = norm(G, 'fro')/sqrt(m*n);
  if params.verbose >= 2
    fprintf('AM %-3d %f %f\n', k, f, normG);
  end

  if normG < params.am_tol
    break;
  end

  % Line search w/ Armijo rule
  alpha = 1;
  for s=1:params.am_ls_max_iters
    Ktemp = K + alpha*Kt;
    max_eig = max(real(eig(A + B*Ktemp)));

    if max_eig < 0
      Ltemp = lyap(A + B*Ktemp, W);
      ftemp = trace(Ltemp*(Q + K'*R*Ktemp)) + ...
              (rho/2)*norm(Ktemp - U, 'fro')^2;
      delta = trace(G'*Kt);
      if ftemp <= f + alpha*params.sigma*delta
        break
      end
    end

    if params.verbose >= 3
      fprintf('AM LS %-3d\n', s);
    end

    alpha = params.beta*alpha;
  end

  K = Ktemp;
  f = ftemp;
  L = Ltemp;
end
