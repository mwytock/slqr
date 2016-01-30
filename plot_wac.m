f1 = h_newton.objval(h_newton.nu < 0);
f2 = h_fista.objval(h_fista.nu < 0);
%f3 = h_admm.objval(h_admm.nu < 0);

min_x = min([f1 f2]);
figure;
semilogy(f1 - min_x);
hold on;
semilogy(f2 - min_x, 'Color', [0 0.5 0]);
%semilogy(f3 - min_x, 'r');

% figure;
% semilogy(h_newton.time(h_newton.nu < 0), f1 - min_x);
% hold on;
% semilogy(h_fista.time(h_fista.nu < 0), f2 - min_x, 'Color', [0 0.5 0]);
% semilogy(h_admm.time(h_admm.nu < 0), f3 - min_x, 'r');