function [A,B,W,Q,R,K_mask] = problem_springs(n)

I = eye(n,n);
Z = zeros(n,n);
T = toeplitz([2 -1 zeros(1,n-2)]);
A = [Z I; -T Z];
B1 = [Z; I];
W = B1*B1';
B = [Z; I];
Q = eye(2*n);
R = 10*I;

K_mask = [I I];