G = @(x,y)1./(1+100*(x.^2-y.^2).^2);
K = chebfun2(G);
F = chebop(@(u) fred(G, u));

l = 0.01;
gaussian_01 = chebfun2(@(x,y) exp(-(x-y).^2/(2*l^2)));

n = 50;
Fn = nystrom(F, gaussian_01, n);