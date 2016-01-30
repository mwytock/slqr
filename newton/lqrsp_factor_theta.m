function X = lqrsp_factor_theta(Theta, theta_eps)
[Z,S,T] = svd(Theta);
s = diag(S);
zr = sqrt(conj(Z(1,:))./T(1,:));
r = min(find((cumsum(s)/sum(s)) > 1-theta_eps));
X = Z(:,1:r)*diag(zr(1:r))*sqrt(S(1:r,1:r));
