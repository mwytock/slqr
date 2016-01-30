% Test solutions for solving the Lyapunov equations
% AX + XA^T + Q = 0

n = 10;
randn('seed', 1);
Q = randn(n);
Q = Q + Q';
A = randn(n);
A = A - 1.1*max(real(eig(A)))*eye(n);

[U,D] = eig(A);
invU = inv(U);
X = lyap_fast(U, D, invU, Q);
assert(norm(X - lyap(A,Q), 'fro') < 1e-8);

disp('PASSED');


