function [D1,D2] = newton_cd(U,S,V,B,W,Q,R,Lambda,K,L,P,G,params)
% Find a low-rank approximation for Theta such that Theta = XX^T
Theta = 1./(S*ones(size(S)) + ones(size(S))*S);
X = lqrsp_factor_theta(Theta, params.cd_theta_eps);

%X = X(:,1:min(10,size(X,2)));
if params.verbose >= 2
  fprintf('Rank of Theta: %d\n', size(X,2));
end

[is,js] = find(K~=0|abs(G)>Lambda);
active = [is js]-1;

if params.verbose >= 2
  fprintf('Active set size: %d\n', size(active,1));
end

if isfield(params, 'save_running')
  save('running.mat','B','R','K','L','P','U','V','X','Lambda','active', ...
       'params');
end

[D1,D2] = newton_cd_loop(B,R,K,L,P,U,V,X,Lambda,active,params);
