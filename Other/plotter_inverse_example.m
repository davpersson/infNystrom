function plotter_inverse_example(filename)
% This function creates Figure 5

addpath('results')
load(filename) % Load the numerical results

% Set up the figure
figure('units','normalized','outerposition',[0 0 3 2])
tiledlayout(2,3)

% Create Figure 5 (a)
nexttile(1)
plot(chebfun2(solution(:,:,1)))
view(2)
title('$t = 0$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
text(-0.15,1,'(a)','Interpreter','latex','Units','normalized','FontSize',20)
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');

% Create Figure 5 (c)
nexttile(4)
plot(chebfun2(solution(:,:,2)))
view(2)
title('$t = 0.1$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
text(-0.15,1,'(c)','Interpreter','latex','Units','normalized','FontSize',20)
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');

% Create Figure 5 (d)
nexttile(5)
plot(chebfun2(solution(:,:,3)))
view(2)
title('$t = 0.5$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
text(-0.15,1,'(d)','Interpreter','latex','Units','normalized','FontSize',20)
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');

% Create Figure 5 (e)
nexttile(6)
plot(chebfun2(solution(:,:,4)))
view(2)
title('$t = 1$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
text(-0.15,1,'(e)','Interpreter','latex','Units','normalized','FontSize',20)
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');

% Create Figure 5 (b)
nexttile(2,[1 2])
semilogy(rank_list,error_optimal/sum(s_true),'k','LineWidth',3)
hold on
semilogy(rank_list,error_list/sum(s_true),'Color',[0 0.4470 0.7410],'LineWidth',3)
xlabel('Rank $k$','Interpreter','latex')
ylabel('Relative trace error','Interpreter','latex')
text(-0.1,1,'(b)','Interpreter','latex','Units','normalized','FontSize',20)
legend({'Optimal','Approximation'},...
    'Interpreter','latex','location','southwest')
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
grid on
hold off

print(filename,'-dpng') % Save figure
end