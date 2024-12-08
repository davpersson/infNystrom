function plotter_pretty_function_example(filename)
% This function creates Figure 1

addpath('results')
load(filename) % Load the numerical results

% Set up the figure
figure('units','normalized','outerposition',[0 0 3 2])
tiledlayout(2,3)

% Create Figure 1 (a)
nexttile(1)
surf(U_true*S_true*U_true')
view(2)
title('Pretty function','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
text(-0.15,1,'(a)','Interpreter','latex','Units','normalized','FontSize',20)
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');

% Create Figure 1 (c)
nexttile(4)
C = approx_cell{1,1};
U = C{1,1}; U = U{1,1}; S = C{1,2}; S = S{1,1};
surf(U*S*U')
view(2)
title('$\ell = 1$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
text(-0.15,1,'(c)','Interpreter','latex','Units','normalized','FontSize',20)
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');

% Create Figure 1 (d)
nexttile(5)
C = approx_cell{1,2};
U = C{1,1}; U = U{1,1}; S = C{1,2}; S = S{1,1};
surf(U*S*U')
view(2)
title('$\ell = 0.1$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
text(-0.15,1,'(d)','Interpreter','latex','Units','normalized','FontSize',20)
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');

% Create Figure 1 (e)
nexttile(6)
C = approx_cell{1,3};
U = C{1,1}; U = U{1,1}; S = C{1,2}; S = S{1,1};
surf(U*S*U')
view(2)
title('$\ell = 0.01$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
text(-0.15,1,'(e)','Interpreter','latex','Units','normalized','FontSize',20)
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');

% Create Figure 1 (b)
nexttile(2,[1 2])
semilogy(rank_list,error_optimal/norm(S_true,'fro'),'k','LineWidth',3)
hold on
semilogy(rank_list,error_list(1,:)/norm(S_true,'fro'),'Color',[0 0.4470 0.7410],'LineWidth',3)
semilogy(rank_list,error_list(2,:)/norm(S_true,'fro'),'Color',[0.6350 0.0780 0.1840],'LineWidth',3)
semilogy(rank_list,error_list(3,:)/norm(S_true,'fro'),'Color',[0.16470588235294117 0.5529411764705883 0.403921568627451],'LineWidth',3)
xlabel('Rank $k$','Interpreter','latex')
ylabel('Relative HS-norm error','Interpreter','latex')
text(-0.1,1,'(b)','Interpreter','latex','Units','normalized','FontSize',20)
legend({'Optimal','$\ell = 1$','$\ell = 0.1$','$\ell = 0.01$'},...
    'Interpreter','latex','location','southwest')
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
grid on
hold off

% Save figure
print(filename,'-dpng')



end