function plotter(filename,plot_type,title_for_plot)

if plot_type == 1
    plotter1(filename,title_for_plot)
elseif plot_type == 2
    plotter2(filename,title_for_plot)
end
end

function plotter1(filename,title_for_plot)
addpath('results')
load(filename)

figure('units','normalized','outerposition',[0 0 0.45 0.6])
surf(U_true*diag(g(diag(S_true)))*U_true')
view(2)
title(title_for_plot,'Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
print(append(filename,'_exact'),'-depsc')

figure('units','normalized','outerposition',[0 0 0.45 0.6])
C = approx_cell{1,1};
U = C{1,1}; U = U{1,1}; S = C{1,2}; S = S{1,1};
surf(U*diag(g(diag(S)))*U')
view(2)
title('$\ell = 1$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
print(append(filename,'_l=1'),'-depsc')

figure('units','normalized','outerposition',[0 0 0.45 0.6])
C = approx_cell{1,2};
U = C{1,1}; U = U{1,1}; S = C{1,2}; S = S{1,1};
surf(U*diag(g(diag(S)))*U')
view(2)
title('$\ell = 0.1$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
print(append(filename,'_l=01'),'-depsc')

figure('units','normalized','outerposition',[0 0 0.45 0.6])
C = approx_cell{1,3};
U = C{1,1}; U = U{1,1}; S = C{1,2}; S = S{1,1};
surf(U*diag(g(diag(S)))*U')
view(2)
title('$\ell = 0.01$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
print(append(filename,'_l=001'),'-depsc')

figure('units','normalized','outerposition',[0 0 0.9 0.45])
semilogy(rank_list,error_optimal/norm(S_true,'fro'),'k','LineWidth',3)
hold on
semilogy(rank_list,error_list(1,:)/norm(S_true,'fro'),'Color',[0 0.4470 0.7410],'LineWidth',3)
semilogy(rank_list,error_list(2,:)/norm(S_true,'fro'),'Color',[0.6350 0.0780 0.1840],'LineWidth',3)
semilogy(rank_list,error_list(3,:)/norm(S_true,'fro'),'Color',[0.16470588235294117 0.5529411764705883 0.403921568627451],'LineWidth',3)
xlabel('Rank $k$','Interpreter','latex')
ylabel('Approximation error','Interpreter','latex')
legend({'Optimal','$\ell = 1$','$\ell = 0.1$','$\ell = 0.01$'},...
    'Interpreter','latex','location','southwest')
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
grid on
hold off
print(append(filename,'_error'),'-depsc')



end


function plotter2(filename)
addpath('results')
load(filename)

figure('units','normalized','outerposition',[0 0 0.45 0.6])
surf(U_true*diag(g(diag(S_true)))*U_true')
view(2)
title(title_for_plot,'Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
print(append(filename,'_exact'),'-depsc')

figure('units','normalized','outerposition',[0 0 0.45 0.6])
C = approx_cell{1,1};
U = C{1,1}; U = U{1,1}; S = C{1,2}; S = S{1,1};
surf(U*diag(g(diag(S)))*U')
view(2)
title('$\ell = 1$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
print(append(filename,'_l=1'),'-depsc')

figure('units','normalized','outerposition',[0 0 0.45 0.6])
C = approx_cell{1,2};
U = C{1,1}; U = U{1,1}; S = C{1,2}; S = S{1,1};
surf(U*diag(g(diag(S)))*U')
view(2)
title('$\ell = 0.1$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
print(append(filename,'_l=01'),'-depsc')

figure('units','normalized','outerposition',[0 0 0.45 0.6])
C = approx_cell{1,3};
U = C{1,1}; U = U{1,1}; S = C{1,2}; S = S{1,1};
surf(U*diag(g(diag(S)))*U')
view(2)
title('$\ell = 0.01$','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
print(append(filename,'_l=001'),'-depsc')

figure('units','normalized','outerposition',[0 0 0.9 0.45])
loglog(rank_list,error_optimal/norm(S_true,'fro'),'k','LineWidth',3)
hold on
loglog(rank_list,error_list(1,:)/norm(S_true,'fro'),'Color',[0 0.4470 0.7410],'LineWidth',3)
loglog(rank_list,error_list(2,:)/norm(S_true,'fro'),'Color',[0.6350 0.0780 0.1840],'LineWidth',3)
loglog(rank_list,error_list(3,:)/norm(S_true,'fro'),'Color',[0.16470588235294117 0.5529411764705883 0.403921568627451],'LineWidth',3)
xlabel('Rank $k$','Interpreter','latex')
ylabel('Approximation error','Interpreter','latex')
legend({'Optimal','$\ell = 1$','$\ell = 0.1$','$\ell = 0.01$'},...
    'Interpreter','latex','location','southwest')
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
grid on
hold off
print(append(filename,'_error'),'-depsc')


end