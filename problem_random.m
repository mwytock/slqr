function [A,B,W,Q,R] = problem_random(n)

randn('state', 1);
A = randn(n)/sqrt(n);
B = randn(n)/sqrt(n);
W = randn(n)/sqrt(n);
Q = randn(n)/sqrt(n);
R = randn(n)/sqrt(n);

W = W'*W;
Q = Q'*Q;
R = R'*R;
