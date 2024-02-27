function plotter(filename)

addpath('results')
load(filename)

figure('units','normalized','outerposition',[0 0 1 1])
subplot(3,2,1)
surf(U_true*diag(g(diag(S_true)))*U_true')
view(2)
title('Kernel','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');

subplot(3,2,2)
C = approx_cell{1,1};
U = C{1,1}; U = U{1,1}; S = C{1,2}; S = S{1,1};
surf(U*diag(g(diag(S)))*U')
view(2)
title('Learned kernel ($\ell = 1$)','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');

subplot(3,2,3)
C = approx_cell{1,2};
U = C{1,1}; U = U{1,1}; S = C{1,2}; S = S{1,1};
surf(U*diag(g(diag(S)))*U')
view(2)
title('Learned kernel ($\ell = 0.1$)','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');

subplot(3,2,4)
C = approx_cell{1,3};
U = C{1,1}; U = U{1,1}; S = C{1,2}; S = S{1,1};
surf(U*diag(g(diag(S)))*U')
view(2)
title('Learned kernel ($\ell = 0.01$)','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');

subplot(3,2,[5 6])
semilogy(rank_list,error_optimal/norm(S_true,'fro'),'k-*','LineWidth',3)
hold on
semilogy(rank_list,error_list(1,:)/norm(S_true,'fro'),'b-*','LineWidth',3)
semilogy(rank_list,error_list(2,:)/norm(S_true,'fro'),'r-*','LineWidth',3)
semilogy(rank_list,error_list(3,:)/norm(S_true,'fro'),'g-*','LineWidth',3)

xlabel('Rank','Interpreter','latex')
ylabel('Relative $\|\cdot\|_{HS}$ error','Interpreter','latex')
legend({'Optimal','$\ell = 1$','$\ell = 0.1$','$\ell = 0.01$'},...
    'Interpreter','latex','location','northeastoutside')
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
hold off

print(filename,'-depsc')


end