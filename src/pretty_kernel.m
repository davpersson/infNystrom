G = @(x,y)1./(1+100*(x.^2-y.^2).^2);
K = chebfun2(G);
F = chebop(@(u) fred(G, u));

l = 0.01;
gaussian_01 = chebfun2(@(x,y) exp(-(x-y).^2/(2*l^2)));
%%
n = 100;

Fn = nystrom_psd(K, gaussian_01, n);
fprintf("MT paper: %.2e\n", norm(K-Fn))

Fn = nystrom_psd2(K, gaussian_01, n);
fprintf("Truncated SVD: %.2e\n", norm(K-Fn))
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

