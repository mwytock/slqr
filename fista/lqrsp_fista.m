function [K,h] = lqrsp_fista(A,B,W,Q,R,Lambda,K0,params)
[m,n] = size(K0);
K = K0;
Kp = K0;
tnu = 1;

tic;
for t=1:params.max_iters
  if tnu == 1
    [At,U,S,V,nu] = lqrsp_eig_deflate(A, B*K, params.max_eig_limit);
  else
    At = A;
    nu = 0;
  end

  % Calculate the Lyapunov equations and the gradient
  if params.acceleration
    Ky = K + (tnu-2)/(tnu+1)*(K-Kp);
  else
    Ky = K;
  end

  % L = lyap_fast(U, S, V, W);
  % P = lyap_fast(V.', S, U.', Q + Ky'*R*Ky);
  L = lyap(At + B*Ky, W);
  P = lyap((At + B*Ky)', Q + Ky'*R*Ky);
  J = trace(P*W);

  normK = sum(sum(Lambda.*abs(K)));
  h.nnz(t) = nnz(K);
  h.nu(t) = nu;
  h.objval(t) = trace(lyap((At + B*K)', Q + K'*R*K)*W) + normK;
  h.norm(t) = normK;
  h.time(t) = toc;

  if params.verbose >= 1
    fprintf('FISTA %-3d %5d %e %e\n', t, nnz(K), h.objval(t), nu);
  end

  if h.time(t) > params.max_time 
    if params.verbose >= 1
      fprintf('Time limit exceeded\n');
    end
    break
  end

  % Use line search to find step size
  G = 2*(R*K + B'*P)*L;
  D = -G;
  alpha = 1;
  for s=1:params.max_ls_iters
    Ktemp = st(Ky + alpha*D, alpha*Lambda);
    [Utemp,Stemp] = eig(At + B*Ktemp);
    max_eig = max(real(diag(Stemp)));

    if max_eig < 0
      Gtemp = (Ky - Ktemp)/alpha;
      Vtemp = inv(Utemp);
      %Ptemp = lyap_fast(Vtemp.', Stemp, Utemp.', Q + Ktemp'*R*Ktemp);
      Ptemp = lyap((At + B*Ktemp)', Q + Ktemp'*R*Ktemp);
      Jtemp = trace(Ptemp*W);

      if Jtemp <= J + alpha*(trace(D'*Gtemp) + 0.5*trace(Gtemp'*Gtemp))
        break
      end
    end

    if params.verbose >= 2
      fprintf('LS %f %f', alpha, max_eig);
      if max_eig < 0
        fprintf(' %f', Jtemp);
      end
      fprintf('\n');
    end

    alpha = params.beta*alpha;
  end


  Kp = K;
  if s == params.max_ls_iters
    if params.verbose >= 1
      fprintf('Exceeded max line searches\n');
    end
    break
  else
    K = Ktemp;
    % if norm(Kp-K,'fro')/sqrt(m*n) < params.tol
    %   if params.verbose >= 1
    %     fprintf('Converged after %d iterations.\n', t);
    %   end
    %   break
    % end
  end

  % Problem changed, reset acceleration
  if nu > 0
    Kp = K;
    tnu = 1;
  else
    tnu = tnu + 1;
  end


end
