function plotter_gaussian_process_example(filename)
% This function creates Figure 4

addpath('results')
load(filename) % Load the numerical results

% Set up the figure
figure('units','normalized','outerposition',[0 0 3 2])
tiledlayout(2,2)

% Generate random Gaussian coefficients for the Karhunen-Loeve expansion
vecOmega2 = randn(size(U,1),1);
vecOmega1 = vecOmega2;

% Plot the exact gaussian process
nexttile(3)
coeffs = reshape(U_true_cts*diag(sqrt(diag(S_true_cts)))*U_true_cts'*vecOmega1,[sqrtN sqrtN]);
surf(legpoly(0:sqrtN-1)*coeffs*legpoly(0:sqrtN-1)')
% surf(U_true_cts*diag(sqrt(diag(S_true_cts)))*...
%     Omega*...
%     diag(sqrt(diag(S_true_cts)))*U_true_cts')
view(2)
title('Exact Gaussian Process','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
%text(-0.15,1,'(a)','Interpreter','latex','Units','normalized','FontSize',20)
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');

% Plot the approximate Gaussian process
nexttile(4)
coeffs = reshape(U*diag(sqrt(diag(S)))*U'*vecOmega2,[sqrtN sqrtN]);
surf(legpoly(0:sqrtN-1)*coeffs*legpoly(0:sqrtN-1)')
view(2)
title('Approximated Gaussian Process','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
%text(-0.15,1,'(a)','Interpreter','latex','Units','normalized','FontSize',20)
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');

% Plot the error plot
nexttile(1,[1 2])
semilogy(rank_list,error_optimal/trace(S_true),'k','LineWidth',3)
hold on
semilogy(rank_list,error_list/trace(S_true),'Color',[0 0.4470 0.7410],'LineWidth',3)
xlabel('Rank $k$','Interpreter','latex')
ylabel('Relative trace error','Interpreter','latex')
text(-0.1,1,'(b)','Interpreter','latex','Units','normalized','FontSize',20)
legend({'Optimal','Approximation'},...
    'Interpreter','latex','location','southwest')
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
grid on
hold off

% Save figure
print(filename,'-dpng')

end