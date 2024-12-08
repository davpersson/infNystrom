function [U,S] = nystrom(Afun,UK,SK,r)
% INPUTS:
%   - Afun: Function handle to compute A*x, where A is an integral
%   operator.
%   - UK: Eigenfunctions of the Gaussian process Kernel
%   - SK: Eigenvalues of the Gaussian process Kernel
%   - r: rank of low-rank approximation
% OUTPUTS:
%   - U: Eigenfunctions of Nystrom approximation
%   - S: Eigenvalues of Nystrom approximations

% Use Algorithm 2.1 in https://epubs.siam.org/doi/pdf/10.1137/21M1466244
% Or 16 in https://doi.org/10.1017/S0962492920000021

%Generate Guassian process
Omega = UK*diag(sqrt(diag(SK)))*randn(size(SK,1),r);
Omega = orth(Omega);

% Compute sketch and core matrix
Y = Afun(Omega);
nu = eps*norm(Y,'fro');
Y_nu = Y + nu*Omega;

try
    x = Omega'*Y_nu;
    % project onto symmetric matrix
    x = 0.5*(x+x');
    C = chol(x);
catch
    print("e")
end

% Compute B so that Nystrom approximation = B*B'
B = (C'\Y_nu')';

%Obtain eigenvalue decomposition of Nystrom approximation
[U,S,~] = svd(B);
S = max(0,S.^2-nu*eye(size(S)));

end