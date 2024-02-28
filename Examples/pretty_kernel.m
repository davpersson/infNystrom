G = @(x,y)1./(1+100*(x.^2-y.^2).^2);
K = chebfun2(G);
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
cfg = struct('axis_option',"title={Exact kernel}");
surf3tikz(gcf, 'fig/pretty_kernel',cfg)

%% l = 1
n = 40;
Fn = nystrom_psd(K, gaussian_1, n);
fprintf("MT paper: %.2e\n", norm(K-Fn))
surf(Fn)
view(0,90)
axis square
colorbar
clim([0,1])
colormap parula
xlabel("$x$",Interpreter="latex")
ylabel("$y$",Interpreter="latex")
cfg = struct('axis_option',"title={$\ell=1$}");
surf3tikz(gcf, 'fig/pretty_kernel_1', cfg)

%% l = 0.1
n = 40;
Fn = nystrom_psd(K, gaussian_01, n);
fprintf("MT paper: %.2e\n", norm(K-Fn))
surf(Fn)
view(0,90)
axis square
colorbar
clim([0,1])
colormap parula
xlabel("$x$",Interpreter="latex")
ylabel("$y$",Interpreter="latex")
cfg = struct('axis_option',"title={$\ell=0.1$}");
surf3tikz(gcf, 'fig/pretty_kernel_01', cfg)

%% l = 0.01
n = 40;
Fn = nystrom_psd(K, gaussian_001, n);
fprintf("MT paper: %.2e\n", norm(K-Fn))
surf(Fn)
view(0,90)
axis square
colorbar
clim([0,1])
colormap parula
xlabel("$x$",Interpreter="latex")
ylabel("$y$",Interpreter="latex")
cfg = struct('axis_option',"title={$\ell=0.01$}");
surf3tikz(gcf, 'fig/pretty_kernel_001', cfg)

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
writematrix([N,Optimal,Error_1,Error_01,Error_001],"fig/error_pretty.csv")
