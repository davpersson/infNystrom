function Fn = nystrom_psd(A, K, n)
% NYSTROM 
%
% NYSTROM(A, K, n) returns a rank n approximation to a self-adjoint linear chebop
% operator using n samples from a Gaussian process with covariance kernel K

% Use Algorithm 2.1 in https://epubs.siam.org/doi/pdf/10.1137/21M1466244
% Or 16 in https://doi.org/10.1017/S0962492920000021

% Sample the Gaussian process
Omega = gaussianprocess(K, n);
Omega = orth(Omega);

% Sample A at X and Y
Y = A*Omega;
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
B = (C'\Y_nu')';

[U,S,~] = svd(B);
S = max(0,S.^2-nu*eye(size(S)));

Fn = U*S*U';
end