function x = compute_hadamard(FD,G,MD,N,i,j,transposeG)
x = 0;
for k=1:size(FD,3)
  if transposeG 
    Gk = G(j,:,k)';
  else
    Gk = G(:,j,k);
  end
  x = x + FD(i,:,k)*Gk + MD(i,:,k)*N(:,j,k);
end
x = real(x);