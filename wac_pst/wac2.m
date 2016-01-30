%svm_mgen;
%save data50m_zico.mat;
load data50m_zico.mat;

A = a_mat;
B = b_vr;
n = size(mac_state,1);

[U,S] = eig(A);
Theta = 1./(S*ones(size(S)) + ones(size(S))*S);
semilogy(svd(Theta));

s0 = svd(Theta);
cumsum(s0)/sum(s0)

th = find(mac_state(:,2) == 1);
thd = find(mac_state(:,2) == 2);
Q = zeros(n,n);
Q(th,th) = 500*((1+0.1)*eye(length(th)) - ones(length(th))/length(th));
Q(thd,thd) = 5000*eye(length(thd));
R = eye(size(B,2));

K_local = sparse(mac_state(:,3),1:size(mac_state,1),ones(size(mac_state,1),1));

K = -lqr(A,B,Q,R);
clf; hold on;
plot(eig(A), 'rx');
plot(eig(A+B*K),'bx')

[U,S] = eig(A+B*K);
Theta = 1./(S*ones(size(S)) + ones(size(S))*S);
semilogy(svd(Theta));

s0 = svd(Theta);
cumsum(s0)/sum(s0)