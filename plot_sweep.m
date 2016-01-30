

figure;
semilogx(lambdas, sweep.nu, 'LineWidth', 1);
set(gca, 'xdir', 'reverse');
prepare_figure('sweep_nu.pdf', [8 6], '\lambda', '\nu')

stable = min(find(sweep.nu <0));
J = sweep.objval - sweep.norm;
Jstar = min(J);

figure;
semilogx(lambdas, (J-Jstar)/Jstar, 'LineWidth', 1);
set(gca, 'xdir', 'reverse');
a = axis();
hold on;
plot([lambdas(stable) lambdas(stable)], a(3:4), 'r', 'LineWidth', 1);
prepare_figure('sweep_J.pdf', [8 6], '\lambda', '(J-J^*)/J^*')

figure;
semilogx(lambdas, sweep.nnz, 'LineWidth', 1);
set(gca, 'xdir', 'reverse');
a = axis();
hold on;
plot([lambdas(stable) lambdas(stable)], a(3:4), 'r', 'LineWidth', 1);
prepare_figure('sweep_nnz.pdf', [8 6], '\lambda', '||K||_0')