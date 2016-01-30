function run_wac2(name, lqrsp, params)

lambdas0 = [logspace(-1, log10(3), 3)];
cases = {'case39'};

for i=1:length(cases);
  lambdas = lambdas0(i,:);
  [A,B,W,Q,R] = problem_wac2();
  %K0 = zeros(size(B'));
  K0 = -lqr(A,B,Q,R);
  for j=2 %1:length(lambdas)
    Lambda = lambdas(j)*ones(size(B'));
    [K{j},h{j}] = lqrsp(A,B,W,Q,R,Lambda,K0,params);
  end
  save(['wac2_' name '_' cases{i} '.mat'], 'lambdas', 'K', 'h');
end
