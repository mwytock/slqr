function [K,h] = run_wac(name, lqrsp, params, offdiag)

cases = {'pst_nedorfler', 'pst_14', 'pst_16', 'pst_ne', 'pst_48', 'pst_50', 'pst_zico50'};
lambdas0 = [logspace(2, log10(0.1353), 3)
            14.1747 0.1485 0.0278
            100 3.8535 0.1024
            100 3.8535 0.1024
            8.1113 0.5462 0.1789
            100 5.5908 0.0705
            100 5.5908 0.0705];
times = [100 100 100 100 1000 1000 3000];

Kpolish = {};
for i=1:length(cases)
  lambdas = lambdas0(i,:);
  params.max_time = times(i);  
  [A,B,W,Q,R,K_mask] = problem_wac(cases{i});
  K0 = -lqr(A,B,Q,R);
  
  for j=1:length(lambdas)

    fprintf('%s lambda=%f numel=%d mask=%d\n', cases{i}, lambdas(j), numel(K0), nnz(K_mask));

    Lambda = (1-K_mask)*lambdas(j);
    if offdiag 
      Lambda = lambdas(j)*(1-K_mask);
    else
      Lambda = lambdas(j)*ones(size(K_mask));
    end

    % Start from FISTA
    clear fparams;
    fparams.max_iters = 1;
    fparams.verbose = 1;
    fparams.max_eig_limit = 0;
    fparams.beta = 0.5;
    fparams.sigma = 1e-4;
    fparams.tol = 1e-5;
    fparams.max_ls_iters = 100;
    fparams.allow_unstable = true;
    fparams.acceleration = false;
    fparams.max_time = Inf;
    K1 = lqrsp_fista(A,B,W,Q,R,Lambda,K0,fparams);

    [K{j},h{j}] = lqrsp(A,B,W,Q,R,Lambda,K1,params);
  end

  filename = ['wac_' name '_' cases{i}];
  if offdiag
    filename = [filename '_offdiag'];
  end

  save(filename, 'lambdas', 'K', 'h');
end
