options = struct('method','l1','gamval',logspace(-4,0,40), ...
                 'rho',100,'maxiter',100,'blksize',[1 1]);
solpath = lqrsp(A,B1,B2,Q,R,options);