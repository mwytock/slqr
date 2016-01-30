%% Weighted sum of Frobenius norm for block sparsity
% 
function G = block_slog(V,gam,rho,blksize,sub_mat_size)

    G = zeros(size(V));
    mm = blksize(1);
    nn = blksize(2);
    eps = 1.e-3;
    
    p = sub_mat_size(1);
    q = sub_mat_size(2);

    for i = 1:p
        for j = 1:q
            
            Vij = V( mm * (i-1) + 1:mm*i, nn * (j-1) + 1:nn*j);
            nVij = norm(Vij,'fro');
            
            rplus  = ( 1/(2*nVij) ) * ( nVij - eps + sqrt( (nVij + eps)^2 - 4 * (gam/rho) ) );
            rminus = ( 1/(2*nVij) ) * ( nVij - eps - sqrt( (nVij + eps)^2 - 4 * (gam/rho) ) );
            
            c = (gam/rho)/eps;
            
            if nVij > c
                G( mm * (i-1) + 1:mm*i, nn * (j-1) + 1:nn*j ) = rplus * Vij;
            else
                
                % form G+, G-, G0
                Gplus  = rplus * Vij;
                Gminus = rminus * Vij;
                Gnot   = 0 * Vij;
                
                % compute the objective values for G+, G-, G0
                phiplus  = gam * log( 1 + (1/eps) * norm(Gplus,'fro') ) +...
                    (rho/2) * norm( Gplus - Vij, 'fro' )^2;
                phiminus = gam * log( 1 + (1/eps) * norm(Gminus,'fro') ) +...
                    (rho/2) * norm( Gminus - Vij, 'fro' )^2;
                phinot   = gam * log( 1 + (1/eps) * norm(Gnot,'fro') ) +...
                    (rho/2) * norm( Gnot - Vij, 'fro' )^2;
                
                % compare the objective values to determine the global minimum
                if (phiplus <= phiminus) && (phiplus <= phinot)
                    G( mm * (i-1) + 1:mm*i, nn * (j-1) + 1:nn*j) = Gplus;            
                end
                if (phiminus <= phiplus) && (phiminus <= phinot)
                    G( mm * (i-1) + 1:mm*i, nn * (j-1) + 1:nn*j) = Gminus;            
                end
                if (phinot <= phiminus)  && (phinot <= phiplus)
                    G( mm * (i-1) + 1:mm*i, nn * (j-1) + 1:nn*j) = Gnot;            
                end
            end              
        end
    end

end