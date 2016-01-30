%% Conjugate gradient method to compute newton direction

function D = newton_cg(A,U,S,V,B,W,Q,R,Lambda,K,L,P,G,params)

E = P*B + K'*R;

% initialization
D = zeros(size(K));
Rk = -G.*(Lambda==0);
Pk = Rk;
normG = norm(G.*(Lambda==0), 'fro');

% conjugate gradient scheme
for t=1:params.max_cg_iters
  % compute H evaluated at Pk
  G1 = B*Pk*L;
  G2 = E*Pk;
  if params.use_lyap
    Lt = lyap(A + B*K, G1+G1');
    Pt = lyap((A + B*K)', G2+G2');
  else
    Lt = lyap_fast(U, S, V, G1+G1');
    Pt = lyap_fast(V', conj(S), U', G2+G2');
  end
  H = 2*((R*Pk + B'*Pt)*L + (R*K + B'*P)*Lt).*(Lambda==0);

  traceHP = trace(H'*Pk);
  normR = norm(Rk,'fro')^2;

  if params.verbose >= 2
    fprintf('CG %-3d ||Rk|| = %e, d2 = %e\n', t, norm(Rk,'fro'), traceHP);
  end

  if traceHP < 0
    if t == 1
      D = -G;
    end
    break
  end

  alpha = normR / traceHP;
  D = D + alpha*Pk;
  Rk = Rk - alpha*H;
  beta = norm(Rk,'fro')^2/normR;
  Pk = Rk + beta*Pk;

  if norm(Rk,'fro') < min(0.5, sqrt(normG)) * normG
    break
  end
end
