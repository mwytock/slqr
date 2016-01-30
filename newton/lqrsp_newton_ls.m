function [Ktemp,Jtemp] = lqrsp_newton_ls(A,B,W,Q,R,Lambda,J,G,K,D,params)

normK = sum(sum(Lambda.*abs(K)));
Jtemp = Inf;
Ktemp = K;

if any(isinf(vec(D))) || any(isnan(vec(D))) || all(vec(D)==0)
  return
end

alpha = 1;
for k=1:params.ls_max_iters
  Ktemp = K + alpha*D;
  [Utemp,Stemp] = eig(A + B*Ktemp);
  max_eig = max(real(diag(Stemp)));

  if max_eig < 0
    normKtemp = sum(sum(Lambda.*abs(Ktemp)));
    % TODO(mwytock): Figure out problem here?
    %Ltemp = lyap_fast(Utemp, Stemp, inv(Utemp), W);
    Ltemp = lyap(A + B*Ktemp, W);
    Jtemp = trace(Ltemp * (Q + Ktemp'*R*Ktemp)) + normKtemp;
    delta = trace(G'*D) + normKtemp - normK;
    if Jtemp <= J + alpha*params.sigma*delta
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

if params.verbose >= 2
  fprintf('LS %f %f %f\n', alpha, max_eig, Jtemp);
end
