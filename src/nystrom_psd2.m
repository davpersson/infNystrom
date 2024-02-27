function Fn = nystrom_psd2(A, K, n)
% INPUTS:
%   - Afun: Function handle to compute A*x, where A is an integral
%   operator.
%   - UK: Eigenfunctions of the Gaussian process Kernel
%   - SK: Eigenvalues of the Gaussian process Kernel
%   - r: rank of low-rank approximation
% OUTPUTS:
%   - U: Eigenfunctions of Nystrom approximation
%   - S: Eigenvalues of Nystrom approximations

%Generate Guassian process
Omega = gaussianprocess(K, n);
Omega = orth(Omega);

% Compute sketch and core matrix
Y = A*(Omega);
M = Omega'*Y;

% Compute eigenvalue decomposition of core matrix
[V,D] = svd(M,'econ'); 
D(D < 5e-16*D(1,1)) = 0; % Truncate very small eigenvalues to 0;

% Compute B so that Nystrom approximation = B*B'
B = Y*V*pinv(diag(sqrt(diag(D))));

%Obtain eigenvalue decomposition of Nystrom approximation
[U,S,~] = svd(B);
S = S^2;
Fn = U*S*U';
end