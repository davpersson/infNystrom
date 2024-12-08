function gaussian_process_example()
% This function runs the numerical experiments for Figure 4.

% --- Defining parameters ---
sqrtN = 100; % The square root of the number of basis vectors for discretization
sqrtn = 25; % The square root of the rank of the covariance kernel
l = 0.4; % Length parameter
tol = 1e-15; % Tolerance for singular values
rank_list = 1:1:30; % List of ranks
filename = 'results/gaussian_process';
% ------

% --- Set up experiment ---

% Discretize in one variable
K = @(x,y) exp(-((x-y).^2)/(2*l^2))/(l*sqrt(2*pi)); % Defining the SE kernel
F = chebfun2(K); % Creating a chebfun2
[U_true_cts,S_true_cts,~] = svd(F);
Q_small = legpoly(0:sqrtN-1); % Generate sqrtN number of legendre polynomials
B = Q_small'*F*Q_small; % Discretize kernel in one variable
[U_true,S_true,~] = svd(B,'econ'); % Compute SVD of discretization in one variable
r = rank(S_true,tol); % Compute the numerical rank of B
U_true = U_true(:,1:r); S_true = S_true(1:r,1:r); % Truncate

% Discretize in two variables
% A = kron(B,B); % The discretized version of the kernel in two variables
U_true = kron(U_true,U_true);
S_true = kron(S_true,S_true);
Afun = @(X) U_true*(S_true*(U_true'*X));
V = zeros(sqrtN^2,sqrtn^2); V(1:sqrtn^2,1:sqrtn^2) = eye(sqrtn^2); % Eigenvectors for Gaussian process kernel
Lambda = eye(sqrtn^2); % Eigenvalues for Gaussian process kernel
% ------

% --- Run experiment ---
error_list = zeros(1,length(rank_list)); % Allocate space for errors

% Compute optimal errors
s_true = sort(diag(S_true),'descend'); % List of singular values in descending order
error_optimal = sum(s_true) - cumsum(s_true); % Trace errors
error_optimal = error_optimal(1:rank_list(end));

rank_counter = 0; % Iteration counter
for r = rank_list
    
    rank_counter = rank_counter + 1 %Update counter

    [U,S] = nystrom(Afun,V,Lambda,r); % Run the Nyström approximation. Will save the last eigenvalue decomposition
    error_list(rank_counter) = sum(s_true) - trace(S);
    
end
%------

U_true_cts = U_true;
S_true_cts = S_true;
% --- Save results ---
save(filename,'error_list','rank_list','error_optimal','U_true_cts','S_true_cts','S_true','F','U','S','sqrtN');
% ------

% --- Plot results ---
plotter_gaussian_process_example(filename)
% -------

end