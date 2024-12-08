% Code to reproduce the experiment in Fig. 4 of the paper

% Define integral operator
k = @(x,y)besselj(0,100*(x.*y+y.^2));
%k = @(x,y)airy(-13*(x.^2.*y+y.^2));

% Define number of experiments
nmin = 1;
nmax = 100;

K = chebfun2(k);
F = chebop(@(u) fred(k, u));
F_adj = chebop(@(u) fred(@(x,y)k(y,x), u));

% Compute SVD of the operator
[U,S,V] = svd(K);
S = diag(S);

% Define covariance kernel
lambda_n = (1:500).^(-3);

% Gaussian kernels
l = 1;
gaussian_1 = chebfun2(@(x,y) exp(-(x-y).^2/(2*l^2)));
l = 0.1;
gaussian_01 = chebfun2(@(x,y) exp(-(x-y).^2/(2*l^2)));
l = 0.01;
gaussian_001 = chebfun2(@(x,y) exp(-(x-y).^2/(2*l^2)));

% Loop over the rank
norm_jacobi = zeros(nmax, 1);
norm_gaussian_1 = zeros(nmax, 1);
norm_gaussian_01 = zeros(nmax, 1);
norm_gaussian_001 = zeros(nmax, 1);

for n = nmin:1:nmax      
    % Lambda_n
    Fn = randomsvd(F, F_adj, lambda_n, n);
    norm_jacobi(n) = norm(K-Fn);

    % Gaussian_1
    Fn = randomsvd(F, F_adj, gaussian_1, n);
    norm_gaussian_1(n) = norm(K-Fn);

    % Gaussian_01
    Fn = randomsvd(F, F_adj, gaussian_01, n);
    norm_gaussian_01(n) = norm(K-Fn);

    % Gaussian_001
    Fn = randomsvd(F, F_adj, gaussian_001, n);
    norm_gaussian_001(n) = norm(K-Fn);

    sprintf("n = %d, Lambda_n = %e, Gaussian_1 = %e, Gaussian_01 = %e, Gaussian_001 = %e", n, norm_jacobi(n), norm_gaussian_1(n), norm_gaussian_01(n), norm_gaussian_001(n))
end

save("exp2.mat", "norm_jacobi", "norm_gaussian_1", "norm_gaussian_01","norm_gaussian_001")

best_bound = zeros(nmax,1);
for n = 1:100
    best_bound(n) = sum(S(n+1:end).^2)^(0.5);
end
best_bound = best_bound + 1e-16;
csvwrite("error_best_svd.csv",[(1:100)',best_bound'])

