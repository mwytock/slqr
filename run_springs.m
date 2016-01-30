function run_springs(name, lqrsp, params, offdiag)

ns = [25 50 250 500];
lambdas0 = [10 1 1e-1
            10 1 1e-1
            1 1e-1 1e-2
            1 1e-1 1e-2];
times = [10 100 3000 6000];

Kpolish = {};
for i=1:length(ns)
  lambdas = lambdas0(i,:);
  params.max_time = times(i);
  [A,B,W,Q,R,K_mask] = problem_springs(ns(i));
  K0 = -lqr(A,B,Q,R);

  for j=1:length(lambdas)
    fprintf('n=%d\tlambda=%f\n', ns(i), lambdas(j));

    if offdiag 
      Lambda = lambdas(j)*(1-K_mask);
    else
      Lambda = lambdas(j)*ones(size(K_mask));
    end

    % Start from FISTA
    clear fparams;
    fparams.max_iters = 1;
    fparams.verbose = 1;
    fparams.max_eig_limit = -1e-1;
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
  
  filename = ['springs_' name '_' num2str(ns(i))];
  if offdiag
    filename = [filename '_offdiag'];
  end

  save(filename, 'lambdas', 'K', 'h');
end
