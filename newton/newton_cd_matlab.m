function D = newton_cd_matlab(U,S,V,B,W,Q,R,Lambda,K,L,P,G,params)
% Find the regularized Newton direction using coordinate descent
[m,n] = size(K);

Theta = 1./(S*ones(size(S)) + ones(size(S))*S);
X = lqrsp_factor_theta(Theta, params.cd_theta_eps);

% Coordinate descent using an explicit form of the Hessian
[is,js] = find(K~=0|abs(G)>Lambda);

E = P*B + K'*R;
VB = V*B;
VL = V*L;
VE = V*E;
UE = U.'*E;
LE = L*E;

% Update the products
r = size(X,2);
F = zeros([m,m,r]);
G = zeros([n,n,r]);
M = zeros([m,n,r]);
N = zeros([m,n,r]);
for i=1:r
  F(:,:,i) = UE.'*diag(X(:,i))*VB;
  G(:,:,i) = VL.'*diag(X(:,i))*U.';
  M(:,:,i) = UE.'*diag(X(:,i))*VL;
  N(:,:,i) = VB.'*diag(X(:,i))*U.';
end

D   = zeros(m,n);
RD  = zeros(m,n);

FD  = zeros(m,n,r);
FTD = zeros(m,n,r);
MDT = zeros(m,m,r);
NDT = zeros(m,m,r);

for t=1:params.cd_max_iters
  Dold = D;
  for idx=1:length(is)
    i = is(idx); j = js(idx);

    a = 2*(R(i,i)*L(j,j) - compute_hadamard_quadratic(F,G,M,N,i,j));
    b = 2*(LE(j,i) + RD(i,:)*L(:,j) - ...
        compute_hadamard(FD, G,MDT,N,i,j,false) - ...
        compute_hadamard(FTD,G,NDT,M,i,j,true));
    c = K(i,j) + D(i,j);

    if params.verbose >= 3
      fprintf('i=%d j=%d a=%f b=%f c=%f\n', i, j, a, b, c);
    end

    mu = -c + st(c - b/a, Lambda(i,j)/a);
    D(i,j) = D(i,j) + mu;

    % Update products
    RD(:,j) = RD(:,j) + mu*R(:,i);
    for k=1:r
      FD(:,j,k)  = FD(:,j,k)  + mu*F(:,i,k);
      FTD(:,j,k) = FTD(:,j,k) + mu*F(i,:,k)';
      MDT(:,i,k) = MDT(:,i,k) + mu*M(:,j,k);
      NDT(:,i,k) = NDT(:,i,k) + mu*N(:,j,k);
    end
  end

  normD = norm(D)/sqrt(m*n);
  diffD = norm(D - Dold)/sqrt(m*n);
  if params.verbose >= 2
    fprintf('CD %-3d %e\n', t, diffD);
  end

  if diffD < params.cd_tol*normD
    break
  end
end
