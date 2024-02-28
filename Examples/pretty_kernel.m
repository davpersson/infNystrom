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
set(gca,'TickLabelInterpreter','latex')
%cfg = struct('colorbar_log',true);
surf3tikz(gcf, 'fig/pretty_kernel')

%% l = 1
n = 50;
Fn = nystrom_psd(K, gaussian_1, n);
fprintf("MT paper: %.2e\n", norm(K-Fn))
surf(Fn)
view(0,90)
axis square
colorbar
clim([0,1])
colormap parula
surf3tikz(gcf, 'fig/pretty_kernel_1')

%% l = 0.1
n = 50;
Fn = nystrom_psd(K, gaussian_01, n);
fprintf("MT paper: %.2e\n", norm(K-Fn))
surf(Fn)
view(0,90)
axis square
colorbar
clim([0,1])
colormap parula
surf3tikz(gcf, 'fig/pretty_kernel_01')

%% l = 0.01
n = 50;
Fn = nystrom_psd(K, gaussian_001, n);
fprintf("MT paper: %.2e\n", norm(K-Fn))
surf(Fn)
view(0,90)
axis square
colorbar
clim([0,1])
colormap parula
surf3tikz(gcf, 'fig/pretty_kernel_001')


%%
n = 100;
MT = [];
TSVD = [];
for i = 1:10
    sprintf("%d / 10",i)
    Fn = nystrom_psd(K, gaussian_01, n);
    MT = [MT, norm(K-Fn)];

    Fn = nystrom_psd2(K, gaussian_01, n);
    TSVD = [TSVD, norm(K-Fn)];
end
fprintf("MT paper: %.2e, TSVD: %.2e\n", mean(MT), mean(TSVD))

