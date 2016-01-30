function D = newton_cd_explicit(A,B,W,Q,R,lambda,K,L,P,gradK,Z,params)
% Find the regularized Newton direction using coordinate descent

% Compute the penalty term on K by looking at the minimum eigenvalue of H
% TODO(mwytock): Calculate the Hessian efficiently and restrict it just to
% the terms in the active set.
G = @(x) vec(lqrsp_gradient(A,B,W,Q,R,reshape(x,size(K))));
H = numdiff(G, vec(K));

% Coordinate descent using an explicit form of the Hessian
g = vec(gradK);
d = zeros(size(g));
is = find(vec(K)~=0|abs(g)>lambda);
k = length(is);

for t=1:params.cd_max_iters
  d_old = d;
  for k=1:numel(is)
    i = is(k);
    a = H(i,i);
    b = H(i,:)*d + g(i);
    c = K(i) + d(i);
    mu = -c + st(c - b/a, lambda/a);
    d(i) = d(i) + mu;
  end
  
  normd = norm(d)/sqrt(k);
  diffd = norm(d - d_old)/sqrt(k);
  if params.verbose >= 2
    fprintf('CD %-3d %e\n', t, diffd);
  end
  
  if diffd < params.cd_tol*normd
    break
  end
end

D = reshape(d,size(K));