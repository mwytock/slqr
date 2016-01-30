function [A,U,S,V,delta] = lqrsp_eig_deflate(A, BK, limit)
n = size(A,1);

[U,S] = eig(A+BK);
V = inv(U);
max_eig = max(real(diag(S)));
delta = max_eig - limit;
if delta > 0
  A = A - delta*eye(n);
  S = S - delta*eye(n);
end
