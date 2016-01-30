%% G-minimization step
%
% truncation operator for cardinality function
function  G = trun(V,gam,rho)           
          b = sqrt(2 * gam / rho);
          G = V .* ( abs(V) > b );
end

