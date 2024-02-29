G = @(x,y)(1+sqrt(5)*abs(x-y)+(5/3)*(x-y).^2).*exp(-sqrt(5)*abs(x-y));
K = chebfun2(G,[500,500]);
l = 1;
gaussian_1 = chebfun2(@(x,y) exp(-(x-y).^2/(2*l^2)));
l = 0.1;
gaussian_01 = chebfun2(@(x,y) exp(-(x-y).^2/(2*l^2)));
l = 0.01;
gaussian_001 = chebfun2(@(x,y) exp(-(x-y).^2/(2*l^2)));

%% Plot kernel
surf(K)
view(0,90)
axis square
colorbar
clim([0,1])
colormap parula
xlabel("$x$",Interpreter="latex")
ylabel("$y$",Interpreter="latex")
set(gca,'TickLabelInterpreter','latex')
cfg = struct('axis_option',"title={$G_{5/2}$}");
surf3tikz(gcf, 'fig/matern_5_2',cfg)

%% Convergence

% Compute svd
[U,S,~] = svd(K);
S = diag(S);
Optimal = [];
Error_1 = [];
Error_01 = [];
Error_001 = [];
N = [];

for i = 5:5:100
    sprintf("%d",i)
    
    N = [N; i];
    % Compute optimal
    Optimal = [Optimal; norm(S(i+1:end))];
    Fn = nystrom_psd(K, gaussian_1, i);
    Error_1 = [Error_1; norm(K-Fn)];
    Fn = nystrom_psd(K, gaussian_01, i);
    Error_01 = [Error_01; norm(K-Fn)];
    Fn = nystrom_psd(K, gaussian_001, i);
    Error_001 = [Error_001; norm(K-Fn)];
end
writematrix([N,Optimal,Error_1,Error_01,Error_001],"fig/error_G_52.csv")
