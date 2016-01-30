% Block sparsity example

% number of systems
N = 5;

% size of each system
nn = 3; 
mm = 1;

% use cyclic condition to obtain unstable system
a = ones(1,nn);
b = 1.5*(sec(pi/nn))*a;

% state-space representation of each system
Aa = -diag(a) + diag(b(2:nn),-1); 
Aa(1,nn) = -b(1);
Bb1 = diag(b);
Bb2 = zeros(nn,1); 
Bb2(1) = b(1);

% non-symmetric weighted Laplacian matrix 

% adjacency matrix 
Ad = toeplitz([1 0 0 1 0 0 1 0 0 1 0 0 1 0 0]);
for i = 1 : N
    for j = 1 : N
        if i ~= j
            cij = 0.5 * ( i - j );
        else
            cij = 0;
        end
        Ad( nn*(i-1)+1 : nn*i, nn*(j-1)+1 : nn*j) = cij * eye(nn);
    end
end
        
% take the sum of each row
d = sum(Ad,2);

% form the Laplacian matrix
L = Ad - diag(d);

% state-space representation of the interconnected system

A  = kron(eye(N), Aa) - L;
B1 = kron(eye(N), Bb1);
B2 = kron(eye(N), Bb2);
Q  = eye(nn*N);
R  = eye(N);

% compute block sparse feedback gains
options_blkwl1 = struct('method','blkwl1','gamval', ...
    logspace(-1,log10(5),50),'rho',100,'maxiter',1000,'blksize',[1 3], ...
    'reweightedIter',1);

tic
solpath_blkwl1 = lqrsp(A,B1,B2,Q,R,options_blkwl1);
toc

% compute element sparse feedback gains
options_wl1 = struct('method','wl1','gamval', ...
    logspace(-1,log10(5),50),'rho',100,'maxiter',1000,'blksize',[1 1], ...
    'reweightedIter',1);

tic
solpath_wl1 = lqrsp(A,B1,B2,Q,R,options_wl1);
toc

% Computational results

% number of nonzero blocks vs. gamma
figure
semilogx(solpath_blkwl1.gam,solpath_blkwl1.nnz,'o','MarkerSize',10,'LineWidth',2)
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', 18)
xlab = xlabel('\gamma','interpreter', 'tex');
set(xlab, 'FontName', 'cmmi10', 'FontSize', 26)

% H2 performance vs. gamma
figure
semilogx(solpath_blkwl1.gam,solpath_blkwl1.Jopt,'r+','MarkerSize',10,'LineWidth',2)
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', 18)
xlab = xlabel('\gamma','interpreter', 'tex');
set(xlab, 'FontName', 'cmmi10', 'FontSize', 26)

% Compare the sparse feedback gains 

% block sparse feedback gain
idx_blk = 46;
solpath_blkwl1.gam(idx_blk)
solpath_blkwl1.nnz(idx_blk)

figure,
spy(solpath_blkwl1.Fopt(:,:,idx_blk),30)
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', 18)
xlabel('')

% element sparse feedback gain
idx = 33;
solpath_wl1.gam(idx)
solpath_wl1.nnz(idx)
figure,spy(solpath_wl1.Fopt(:,:,idx),30)
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', 18)
xlabel('')

(solpath_blkwl1.Jopt(idx_blk) - solpath_wl1.Jopt(idx))/solpath_wl1.Jopt(idx)