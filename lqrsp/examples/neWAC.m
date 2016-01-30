% Wide Area Control Example
load neData.mat

% data description:
% A, B1, B2 are state space matrices
% Q and R are the state and control performance weights
% uo is the matrix of right eigenvectors for open-loop
% uc is the matrix of right eigenvectors for closed-loop

% ADMM
options = struct('method','wl1','gamval',logspace(-4,0,40), ...
        'rho',100,'maxiter',100,'blksize',[1 1],'reweightedIter',5);
tic
solpath = lqrsp(A,B1,B2,Q,R,options);
toc

% Performance vs. Sparsity
figure
% number of nonzeros of feedback gains vs. gamma
semilogx(solpath.gam,solpath.nnz/675*100,'.','LineWidth',2,'MarkerSize',25);
h = get(gcf,'CurrentAxes');
set(h,'FontSize',18)
xlab = xlabel('\gamma','interpreter', 'tex');
set(xlab,'FontSize', 15);
ylabel('percent','FontName', 'cmmr10','interpreter', 'tex');
set(gca,'XTick',[0.0001,0.001,0.01,0.1,1],'YTick',[0,20,40,60,80]);
title('card(F) / card(F_c)','FontName', 'cmmr10');

% H2 norm vs. gamma
figure
[Fc, P] = lqr(A,B2,Q,R);
Jc = trace(P*(B1*B1'));
semilogx(solpath.gam,(solpath.Jopt - Jc)/Jc*100,...
    'r.','LineWidth',2,'MarkerSize',25)
h = get(gcf,'CurrentAxes');
set(h,'FontSize', 18)
xlab = xlabel('\gamma','interpreter', 'tex');
set(xlab,'FontSize', 15)
ylabel('percent','FontName', 'cmmr10','interpreter', 'tex');
set(gca,'XTick',[0.0001,0.001,0.01,0.1,1],...
'YTick',[0,0.4,0.8,1.2,1.6],'YLim',[0 1.7])
title('(J - J_c) / J_c','FontName', 'cmmr10');

% Compass plots of dominant inter-area modes

% indices of angles and speeds(frequency)
speed_idx = [2,9,17,25,33,41,49,57,65,73]';

% inter-area mode 1(OPEN LOOP): 
% critical values are 
lamc1_1 = 20;
lamc1_2 = 21; 

% identify critical speeds (from speed eigenvector)
crit_set1_1 = [1:9];
crit_set1_2 = [10];

% corresponding states and output
set1_1 = speed_idx(crit_set1_1);
set1_2 = speed_idx(crit_set1_2);

% Compass plot of right eigenvectors of speeds
figure
compass(uo(speed_idx,lamc1_1));
th = findall(gcf,'Type','text');
for i = 1:length(th),
set(th(i),'FontSize',15)
end
hold on;
h1 = compass(uo(set1_1,lamc1_1),'g');
set(h1,'LineWidth',3);
h2 = compass(uo(set1_2,lamc1_1),'r');
set(h2,'LineWidth',3);
title('Mode 1','FontSize',18)

% inter-area mode 1(CLOSED LOOP): 
% critical values are 
lamc1_1c = 25;
lamc1_2c = 26; 
% identify critical speeds (from speed eigenvector)
crit_set1_1c = [1:9];
crit_set1_2c = [10];

% corresponding states and output
set1_1c = speed_idx(crit_set1_1c);
set1_2c = speed_idx(crit_set1_2c);

% Compass plot of right eigenvectors of speeds
figure
compass(-uc(speed_idx,lamc1_1c));
th = findall(gcf,'Type','text');
for i = 1:length(th),
set(th(i),'FontSize',15)
end
hold on;
h3 = compass(-uc(set1_1c,lamc1_1c),'g');
set(h3,'LineWidth',3);
h4 = compass(-uc(set1_2c,lamc1_1c),'r');
set(h4,'LineWidth',3);
title('Mode 1','FontSize',18)

% Signal exchange network
% gamma = 0.0289, nnz = 90
figure
Ko = solpath.Fopt(:,:,25);
Ko(1,1:7) = 0; Ko(2,8:15) = 0;
Ko(3,16:23) = 0; Ko(4,24:31) = 0;
Ko(5,32:39) = 0; Ko(6,40:47) = 0;
Ko(7,48:55) = 0; Ko(8,56:63) = 0;
Ko(9,64:71) = 0;

Kc = solpath.Fopt(:,:,25) - Ko;
spy(Ko,'r'); xlabel('');
hold on
spy(Kc,'b'); xlabel('');
title('\gamma = 0.0289, card(F) = 90','FontSize',18)

% gamma = 1, nnz = 37
figure
Ko = solpath.Fopt(:,:,40);
Ko(1,1:7) = 0; Ko(2,8:15) = 0;
Ko(3,16:23) = 0; Ko(4,24:31) = 0;
Ko(5,32:39) = 0; Ko(6,40:47) = 0;
Ko(7,48:55) = 0; Ko(8,56:63) = 0;
Ko(9,64:71) = 0;

Kc = solpath.Fopt(:,:,40) - Ko;
spy(Ko,'r'); xlabel('');
hold on
spy(Kc,'b'); xlabel('');
title('\gamma = 1, card(F) = 37','FontSize',18)