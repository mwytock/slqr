function G = lqrsp_gradient(A,B,W,Q,R,K)

L = lyap((A + B*K), W);    
P = lyap((A + B*K)', Q + K'*R*K);
G = 2*(R*K + B'*P)*L;
