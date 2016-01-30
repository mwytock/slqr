function J = lqrsp_objective(A,B,W,Q,R,K)

L = lyap(A + B*K, W);
J = trace(L*(Q + K'*R*K));