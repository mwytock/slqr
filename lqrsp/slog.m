%% G-minimization step
%
% proximity operator for sum-of-log function

function G  =  slog(V,gam,rho)

        eps = 1.e-3;
          c = (gam/rho) / eps;
             
         % the solutions of the quadratic equations    
         Gp_plus  = (1/2) * (V - eps + sqrt( (V + eps).^2 - 4 * (gam/rho) ));         
         Gp_minus = (1/2) * (V - eps - sqrt( (V + eps).^2 - 4 * (gam/rho) ));
         Gm_plus  = (1/2) * (V + eps + sqrt( (V - eps).^2 - 4 * (gam/rho) ));
         Gm_minus = (1/2) * (V + eps - sqrt( (V - eps).^2 - 4 * (gam/rho) ));
                  
         % objective function values
         phi_Gp_plus  = phi_log(Gp_plus, V,gam,rho,eps);
         phi_Gp_minus = phi_log(Gp_minus,V,gam,rho,eps);         
         phi_Gm_plus  = phi_log(Gm_plus, V,gam,rho,eps);
         phi_Gm_minus = phi_log(Gm_minus,V,gam,rho,eps);
                 phi0 = phi_log(zeros(size(V)),V,gam,rho,eps);
                 
         % compare the objective function values to determine the global
         % minimum         
         Gp0 = Gp_plus  .* ((phi_Gp_plus <= phi_Gp_minus) & (phi_Gp_plus <= phi0))...
             + Gp_minus .* ((phi_Gp_minus <= phi_Gp_plus) & (phi_Gp_minus <= phi0))...
             +        0 .* ((phi0 <= phi_Gp_plus)         & (phi0 <= phi_Gp_minus));
         Gm0 = Gm_plus  .* ((phi_Gm_plus <= phi_Gm_minus) & (phi_Gm_plus <= phi0))...
             + Gm_minus .* ((phi_Gm_minus <= phi_Gm_plus) & (phi_Gm_minus <= phi0))...
             +        0 .* ((phi0 <= phi_Gm_plus)         & (phi0 <= phi_Gm_minus));
         
         % construct the solution
         G = Gp_plus  .* (V > c)  + Gp0 .* ( (0 <= V) & (V <= c) ) +...
                Gm_minus .* (V < -c) + Gm0 .* ( (-c <= V) & (V < 0) );
end

function phi = phi_log(G,V,gam,rho,eps)

        % Since the solution to the G-minimization problem needs to be
        % real-valued,
        % we set the objective function phi to be infinity for
        % complex-valued elements of G
        
        idxG = double(imag(G)~=0);
         phi = gam * log( 1 + abs(G)/eps ) + (rho/2) * ( G - V ).^2;
         phi = phi .* (1 - idxG) + 1.e+16 * idxG;
         
end

