function plotter_matern_example(filename1,filename2,upper_bound)
% This function creates Figure 2. If upper_bound == 1 it will show the
% plots with the theoretical upper bounds. Otherwise, it wont. 

addpath('results')
load(filename1) % Load the numerical results for Matern-1/2

% Set up the figure
figure('units','normalized','outerposition',[0 0 3 2])
tiledlayout(3,3)

% --- Generate sample paths of Gaussian processes/random fields with the
% covariance kernels in covariance_cells. We will plot these sample paths
% ---

% Unpack covariance cell
C1 = covariance_cell{1,1};
C2 = covariance_cell{1,2};
C3 = covariance_cell{1,3};
V1 = C1{1,1}; V1 = V1{1,1}; Lambda1 = C1{1,2}; Lambda1 = Lambda1{1,1};
V2 = C2{1,1}; V2 = V2{1,1}; Lambda2 = C2{1,2}; Lambda2 = Lambda2{1,1};
V3 = C3{1,1}; V3 = V3{1,1}; Lambda3 = C3{1,2}; Lambda3 = Lambda3{1,1};

% Generate random gaussians
omega1 = randn(size(V1,2),1);
omega3 = randn(size(V3,2),1);
omega2 = randn(size(V2,2),1);

% Karhunen-Loeve expansion
omega1 = V1*diag(sqrt(diag(Lambda1)))*omega1;
omega2 = V2*diag(sqrt(diag(Lambda2)))*omega2;
omega3 = V3*diag(sqrt(diag(Lambda3)))*omega3;

% Normalize the sample paths so that the maximum absolute value is 1
omega1 = omega1/max(abs(omega1));
omega2 = omega2/max(abs(omega2));
omega3 = omega3/max(abs(omega3));
%------

% Create Figure 2 (b)
nexttile(2,[1 2])
semilogy(rank_list,error_optimal/trace(S_true),'k','LineWidth',3)
hold on
semilogy(rank_list,error_list(1,:)/trace(S_true),'Color',[0 0.4470 0.7410],'LineWidth',3)
semilogy(rank_list,error_list(3,:)/trace(S_true),'Color',[0.16470588235294117 0.5529411764705883 0.403921568627451],'LineWidth',3)
semilogy(rank_list,error_list(2,:)/trace(S_true),'Color',[0.6350 0.0780 0.1840],'LineWidth',3)
xlabel('Rank $k$','Interpreter','latex')
ylabel('Relative trace error','Interpreter','latex')
text(-0.08,1.1,'(b)','Interpreter','latex','Units','normalized','FontSize',20)
legend({'Optimal','$\ell = 1$','$\ell = 0.01$','Matern-3/2','','',''},...
    'Interpreter','latex','location','northeast')
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
L = legend;
L.AutoUpdate = 'off';
grid on

% Add theoretical upper bounds
if upper_bound == 1

    ub = upper_bound_generator(U_true,S_true,covariance_cell,rank_list);
    semilogy(rank_list,ub(1,:)/trace(S_true),'--','Color',[0 0.4470 0.7410],'LineWidth',3)
    semilogy(rank_list,ub(3,:)/trace(S_true),'--','Color',[0.16470588235294117 0.5529411764705883 0.403921568627451],'LineWidth',3)
    semilogy(rank_list,ub(2,:)/trace(S_true),'--','Color',[0.6350 0.0780 0.1840],'LineWidth',3)

end

hold off

% Create Figure 2 (a)
nexttile(1)
plot(F)
view(2)
title('Matern-1/2','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
text(-0.15,1.1,'(a)','Interpreter','latex','Units','normalized','FontSize',20)
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');

load(filename2) % Load results for Matern-5/2

% Create Figure 2 (d)
nexttile(5,[1 2])
semilogy(rank_list,error_optimal/trace(S_true),'k','LineWidth',3)
hold on
semilogy(rank_list,error_list(1,:)/trace(S_true),'Color',[0 0.4470 0.7410],'LineWidth',3)
semilogy(rank_list,error_list(3,:)/trace(S_true),'Color',[0.16470588235294117 0.5529411764705883 0.403921568627451],'LineWidth',3)
semilogy(rank_list,error_list(2,:)/trace(S_true),'Color',[0.6350 0.0780 0.1840],'LineWidth',3)
xlabel('Rank $k$','Interpreter','latex')
ylabel('Relative trace error','Interpreter','latex')
text(-0.08,1.1,'(d)','Interpreter','latex','Units','normalized','FontSize',20)
legend({'Optimal','$\ell = 1$','$\ell = 0.01$','Matern-3/2'},...
    'Interpreter','latex','location','northeast')
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
L = legend;
L.AutoUpdate = 'off';
grid on

% Add theoretical upper bounds
if upper_bound == 1

    ub = upper_bound_generator(U_true,S_true,covariance_cell,rank_list);
    semilogy(rank_list,ub(1,:)/trace(S_true),'--','Color',[0 0.4470 0.7410],'LineWidth',3)
    semilogy(rank_list,ub(3,:)/trace(S_true),'--','Color',[0.16470588235294117 0.5529411764705883 0.403921568627451],'LineWidth',3)
    semilogy(rank_list,ub(2,:)/trace(S_true),'--','Color',[0.6350 0.0780 0.1840],'LineWidth',3)

end

hold off

% Create Figure 2 (c)
nexttile(4)
plot(F)
view(2)
title('Matern-5/2','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
text(-0.15,1.1,'(c)','Interpreter','latex','Units','normalized','FontSize',20)
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');

% Create Figure 2 (e)
nexttile(7)
plot(V2*Lambda2*V2')
view(2)
title('Matern-3/2','Interpreter','latex')
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
text(-0.15,1.1,'(e)','Interpreter','latex','Units','normalized','FontSize',20)
colorbar('TickLabelInterpreter', 'latex');
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');

% Create Figure 2 (f)
nexttile(8,[1 2])
plot(omega1,'Color',[0 0.4470 0.7410],'LineWidth',3)
hold on
plot(omega3,'Color',[0.16470588235294117 0.5529411764705883 0.403921568627451],'LineWidth',3)
plot(omega2,'Color',[0.6350 0.0780 0.1840],'LineWidth',3)
xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
text(-0.04,1.1,'(f)','Interpreter','latex','Units','normalized','FontSize',20)
legend({'$\ell = 1$','$\ell = 0.01$','Matern-3/2'},...
    'Interpreter','latex','location','northeast')
set(gca,'FontSize',20)
set(gca,'TickLabelInterpreter','latex');
hold off

% Save figure
print('results/matern_example','-dpng')

end

function ub = upper_bound_generator(U_true,S_true,covariance_cell,rank_list)
% This function creates the upper bounds in Theorem 3.3 (3.6c)

% Initialize a matrix containing the lists of betas
beta_list = zeros(length(covariance_cell),length(rank_list));
delta_list = zeros(length(covariance_cell),length(rank_list));

for i = 1:length(covariance_cell)

    % Unpack covariance_cell
    C = covariance_cell{1,i};
    V = C{1,1}; V = V{1,1}; Lambda = C{1,2}; Lambda = Lambda{1,1};
    sqrtS_true = diag(sqrt(diag(S_true)));
    sqrtLambda = diag(sqrt(diag(Lambda)));

    for j = 1:length(rank_list)

        k = rank_list(j);

        K22 = (U_true(:,k+1:end)'*V) * Lambda * (U_true(:,k+1:end)'*V)';
        K11 = (U_true(:,1:k)'*V) * Lambda * (U_true(:,1:k)'*V)';
        s = svd(K11);
        Q = orth(sqrtLambda * V' * U_true(:,1:k));
        
        delta = norm(sqrtS_true(k+1:end,k+1:end)*(U_true(:,k+1:end)'*V) * sqrtLambda *Q ,'fro')^2;
        beta = (norm(sqrtS_true(k+1:end,k+1:end)*(U_true(:,k+1:end)'*V) * sqrtLambda ,'fro')^2 - delta);
        delta = norm(sqrtS_true(k+1:end,k+1:end)*(U_true(:,k+1:end)'*V) * Lambda * (V'*U_true(:,1:k)) * inv(K11),'fro')^2;


        % If the covariance matrix is numerically low-rank, just assign the
        % value 1e15 to the beta and delta. In reality, beta = delta =
        % inf in these cases
        if max(1./s) > 1e15
            beta = 1e15;
            delta = 1e15;
        else
            beta = beta*max(1./s);
        end
        beta_list(i,j) = beta;
        delta_list(i,j) = delta;


    end


end

% Now when we know the betas and deltas, we can create the upper bounds.
ub = zeros(length(covariance_cell),length(rank_list));

for i = 1:length(covariance_cell)

    ub(i,[1 2]) = trace(S_true)*ones(1,2);

    ub_temp = trace(S_true);

    for j = 3:length(rank_list)

        l = rank_list(j);

        for p = 2:(l-1)

            k = l-p;
            ub_temp = min(ub_temp,trace(S_true(k+1:end,k+1:end)) + delta_list(i,k) + (k/(p-1))*beta_list(i,k));

        end

        ub(i,j) = ub_temp;


    end


end



end