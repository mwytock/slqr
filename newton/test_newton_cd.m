% Test thew Newton-CD method on the 5 mass case

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
Lambda = 1*ones(size(K0));

clear params;
params.max_iters = 1;
params.sigma = 1e-4;
params.beta = 0.5;
params.verbose = 2;
params.proj_orthant = 0;
params.cd_max_iters = 20;
params.cd_tol = 1e-6;
params.cd_theta_eps = 1e-2;

% Test at the initial condition, one iteration
params.cd_max_iters = 1;
params.verbose = 3;
K = K0;
L = lyap((A + B*K), W);
J = trace(L*(Q + K'*R*K)) + sum(sum(Lambda.*abs(K)));
P = lyap((A + B*K)', Q + K'*R*K);
D1 = newton_cd_matlab(A,B,W,Q,R,Lambda,K,L,P,gradK,Z,params);
D2 = newton_cd(A,B,W,Q,R,Lambda,K,L,P,gradK,Z,params);
assert(norm(D1-D2) < 1e-12);

% Initial condition until convergence
params.cd_max_iters = 20;
params.verbose = 2;
D1 = newton_cd_matlab(A,B,W,Q,R,Lambda,K,L,P,gradK,Z,params);
D2 = newton_cd(A,B,W,Q,R,Lambda,K,L,P,gradK,Z,params);
assert(norm(D1-D2) < 1e-12);

% Step 2, one iteration
K = K0 + D1;
params.cd_max_iters = 1;
params.verbose = 3;
D1 = newton_cd_matlab(A,B,W,Q,R,Lambda,K,L,P,gradK,Z,params);
D2 = newton_cd(A,B,W,Q,R,Lambda,K,L,P,gradK,Z,params);
assert(norm(D1-D2) < 1e-12);

% Step 2, until convergence
params.cd_max_iters = 20;
params.verbose = 2;
D1 = newton_cd_matlab(A,B,W,Q,R,Lambda,K,L,P,gradK,Z,params);
D2 = newton_cd(A,B,W,Q,R,Lambda,K,L,P,gradK,Z,params);
assert(norm(D1-D2) < 1e-12);

% Test the whole thing
params.max_iters = 5;
params.verbose = 1;
params.argmin_dir = @newton_cd;
[K,h] = lqrsp_newton(A,B,W,Q,R,Lambda,K0,params);

assert(abs(h.objval(end) - 24.8203) < 1e-4)
assert(abs(K(1,2) - (-0.0205))    < 1e-4);
assert(K(2,1) == 0)

disp('PASSED');
