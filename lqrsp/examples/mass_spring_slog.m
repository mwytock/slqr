% Mass-spring system

% State-space representation of the mass-spring system with N = 50 masses
N = 50;
I = eye(N,N);
Z = zeros(N,N);
T = toeplitz([2 -1 zeros(1,N-2)]);
A = [Z I; -T Z]; 
B1 = [Z; I]; 
B2 = [Z; I]; 
Q = eye(2*N); 
R = 10*I;

% Compute the optimal sparse feedback gains
options = struct('method','slog','gamval',logspace(-4,-1,50),...
        'rho',100,'maxiter',100,'blksize',[1 1]);
tic
solpath = lqrsp(A,B1,B2,Q,R,options);
toc

% Computational results

% number of nonzeros of feedback gains vs. gamma
figure(1)
semilogx(solpath.gam,solpath.nnz,'o','LineWidth',2,'MarkerSize',10);
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', 18)
xlab = xlabel('\gamma','interpreter', 'tex');
set(xlab, 'FontName', 'cmmi10', 'FontSize', 26)

idx = [1,40,50];

% Sparsity patterns vs. gamma
for k = 1:length(idx)
    figure
    spy(solpath.F(:,:,idx(k)),10)
    h = get(gcf,'CurrentAxes');
    set(h, 'FontName', 'cmr10', 'FontSize', 18)        
    xlabel('')
end

% H2 norm vs. gamma
[Fc, P] = lqr(A,B2,Q,R);
Jc = trace(P*(B1*B1'));
figure
semilogx(solpath.gam,(solpath.Jopt - Jc)/Jc*100,...
    'r+','LineWidth',2,'MarkerSize',10)
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', 18, 'yscale', 'log')
xlab = xlabel('\gamma','interpreter', 'tex');
set(xlab, 'FontName', 'cmmi10', 'FontSize', 26)
set(gca,'YTick',[0.1,1,10],'YTickLabel',{'0.1%','1%','10%'}, 'XLim', [0.001 0.1])

% number of nonzero element vs. gamma
figure,
semilogx(solpath.gam,solpath.nnz/(2*N^2)*100,...
    'o','LineWidth',2,'MarkerSize',10)
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', 18, 'yscale', 'lin')
xlab = xlabel('\gamma','interpreter', 'tex');
set(xlab, 'FontName', 'cmmi10', 'FontSize', 26)
set(gca,'YTick',[1,10,20],'YTickLabel',{'1%','10%','20%'}, 'XLim', [0.001 0.1])

% diagonal part of the position and velocity feedback gains
% for different values of gamma

idx = [1 40 50];

figure
hold on
plot(1:N,diag(solpath.Fopt(:,1:N,idx(1))),'o','LineWidth',2,'MarkerSize',10);
plot(1:N,diag(solpath.Fopt(:,1:N,idx(2))),'r+','LineWidth',2,'MarkerSize',10);
plot(1:N,diag(solpath.Fopt(:,1:N,idx(3))),'k*','LineWidth',2,'MarkerSize',10);
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', 18)    
xlab = xlabel('n','interpreter', 'tex');
set(xlab, 'FontName', 'cmmi10', 'FontSize', 26)

figure
hold on
plot(1:N,diag(solpath.Fopt(:,N+1:2*N,idx(1))),'o','LineWidth',2,'MarkerSize',10);
plot(1:N,diag(solpath.Fopt(:,N+1:2*N,idx(2))),'r+','LineWidth',2,'MarkerSize',10);
plot(1:N,diag(solpath.Fopt(:,N+1:2*N,idx(3))),'k*','LineWidth',2,'MarkerSize',10);
h = get(gcf,'CurrentAxes');
set(h, 'FontName', 'cmr10', 'FontSize', 18)    
xlab = xlabel('n','interpreter', 'tex');
set(xlab, 'FontName', 'cmmi10', 'FontSize', 26)

% indices for sparsity vs. performance 
% [30 44 50]
% idx = [30 44 50];
%  solpath.gam(idx) 
%  solpath.nnz(idx)/(2*N^2)*100
%  (solpath.Jopt(idx) - Jc)/Jc * 100