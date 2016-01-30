% Test the derivative/second differential for the LQR problem

N = 5;
I = eye(N,N);
Z = zeros(N,N);
T = toeplitz([2 -1 zeros(1,N-2)]);
A = [Z I; -T Z];
B1 = [Z; I];
W = B1*B1';
B = [Z; I];
Q = eye(2*N);
R = 10*I;

K0 = -lqr(A,B,Q,R);
J = @(K) lqrsp_objective(A,B,W,Q,R,K);

rand('seed', 1);

% At the optimum
K = K0;
L = lyap((A + B*K), W);
P = lyap((A + B*K)', Q + K'*R*K);
gradK = 2*(R*K + B'*P)*L;
gradK2 = numdiff(J, K);
assert(norm(gradK - gradK2, 'fro') < 1e-4);

% At a random point
K(1:4,:) = 0;
L = lyap((A + B*K), W);
P = lyap((A + B*K)', Q + K'*R*K);
gradK = 2*(R*K + B'*P)*L;
gradK2 = numdiff(J, K);
assert(norm(gradK - gradK2, 'fro') < 1e-4);

% The second differential evaluted at the optimum
J = @(x) lqrsp_objective(A,B,W,Q,R,0,reshape(x,size(K0)));

% At the optimum
K = K0;
D = randn(size(K));

L = lyap((A + B*K), W);
P = lyap((A + B*K)', Q + K'*R*K);
Lt = lyap((A + B*K), B*D*L + (B*D*L)');
Pt = lyap((A + B*K)', (K'*R + P*B)*D + D'*(K'*R + P*B)');
H = 2*((R*D + B'*Pt)*L + (R*K + B'*P)*Lt);
[g, H2] = numdiff(J, vec(K));
assert(abs(vec(D)'*H2*vec(D) - trace(H'*D))/norm(H,'fro') < 1e-4);

% At a random point
K(1:4,:) = 0;
D = randn(size(K));

L = lyap((A + B*K), W);
P = lyap((A + B*K)', Q + K'*R*K);
Lt = lyap((A + B*K), B*D*L + (B*D*L)');
Pt = lyap((A + B*K)', (K'*R + P*B)*D + D'*(K'*R + P*B)');
H = 2*((R*D + B'*Pt)*L + (R*K + B'*P)*Lt);
[g, H2] = numdiff(J, vec(K));
assert(abs(vec(D)'*H2*vec(D) - trace(H'*D))/norm(H,'fro') < 1e-4);

% Test computing the Hessian by taking the Jacobian of the gradient
G = @(x) vec(lqrsp_gradient(A,B,W,Q,R,reshape(x,size(K0))));
H3 = numdiff(G, vec(K));
assert(norm(H2-H3, 'fro')/norm(H2,'fro') < 1e-4);

% Test the explicit form of the Lyapunov function
[U,Lambda] = eig(A+B*K);
Theta = 1./(Lambda*ones(size(Lambda)) + ones(size(Lambda))*Lambda);
V = inv(U);
Lt2 = -U*((V*(B*D*L + L*D'*B')*V.') .* Theta)*U.';
Pt2 = -V.'*((U.'*(D'*(B'*P + R*K) + (P*B + K'*R)*D)*U) .* Theta)*V;
assert(norm(Lt - Lt2, 'fro') < 1e-4);
assert(norm(Pt - Pt2, 'fro') < 1e-4);

% Test the explicit form of the second differential using the Theta
% decomposition
[Z,S,T] = svd(Theta);
r = sqrt(conj(Z(1,:))./T(1,:));
X = Z*diag(r)*sqrt(S);
assert(norm(Theta - X*X.') < 1e-8);

% trace \tilde{L}ED and L\tilde{P}BD terms
d = 0;
for k=1:size(X,2)
  Xk = U*diag(X(:,k))*V;
  d = d - trace(Xk.'*D'*E'*Xk*B*D*L) - trace(Xk.'*E*D*Xk*B*D*L);
end
assert(abs(trace(Lt*E*D) - d) < 1e-4);
assert(abs(trace(Pt*B*D*L) - d) < 1e-4);

% The whole thing
d = 2*trace(L*E*D) + trace(L*D'*R*D) + 2*d;
dnum = g'*vec(D) + 0.5*vec(D)'*H2*vec(D);
assert(abs(dnum - d) < 1);

disp('PASSED');
