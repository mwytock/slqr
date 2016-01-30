%% Weighted sum of Frobenius norm for block sparsity
% 
function G = block_trun(V,gam,rho,blksize,sub_mat_size)

    G = zeros(size(V));
    mm = blksize(1);
    nn = blksize(2);
    p = sub_mat_size(1);
    q = sub_mat_size(2);

    for i = 1:p
        for j = 1:q

            Vij = V( mm * (i-1) + 1:mm*i, nn * (j-1) + 1:nn*j);
            nVij = norm(Vij,'fro');

            b = sqrt(2*gam/rho);

            if nVij <= b
                G( mm * (i-1) + 1:mm*i, nn * (j-1) + 1:nn*j) = 0;
            else                 
                G( mm * (i-1) + 1:mm*i, nn * (j-1) + 1:nn*j) = Vij;            
            end              
        end
    end

end