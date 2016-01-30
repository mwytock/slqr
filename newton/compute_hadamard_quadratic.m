function x = compute_hadamard_quadratic(F,G,M,N,i,j)
x = 0;
for k=1:size(F,3)
  x = x + 2*(F(i,i,k)*G(j,j,k) + M(i,j,k)*N(i,j,k));
end
x = real(x);