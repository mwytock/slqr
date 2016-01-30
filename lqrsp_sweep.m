function sweep = lqrsp_sweep(problem, lqrsp, params, lambdas, K0)

[A,B,W,Q,R] = problem();
K = K0;

tic;
for i=1:length(lambdas)
  Lambda = lambdas(i)*ones(size(K));
  fprintf('Iteration %-3d %e\n', i, lambdas(i));
  [K,h] = lqrsp(A,B,W,Q,R,Lambda,K,params);

  sweep.objval(i) = h.objval(end);
  sweep.norm(i) = h.norm(end);
  sweep.nnz(i) = h.nnz(end);
  sweep.nu(i) = h.nu(end);
  sweep.iters(i) = length(h.objval);
  sweep.time(i) = toc;
end
