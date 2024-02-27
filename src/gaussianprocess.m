function f = gaussianprocess(kernel, varargin)
% GAUSSIANPROCESS 
%
% GAUSSIANPROCESS(K, n) returns n samples of the Gaussian process with covariance
% kernel K.
%
% GAUSSIANPROCESS(LAMBDA) samples the Gaussian process using the covariance
% kernel with eigenvectors given by the first l Jacobi polynomials P^(2,2) 
% and eigenvalues Lambda, where l = length(Lambda).

% Get the number of samples
if nargin > 1 && isnumeric(varargin{1})
    n = varargin{1};
else
    n = 1;
end

% Compute the SVD decomposition of the kernel
if isa(kernel, "chebfun2")
    L = chol(kernel, 'lower');
    l = rank(L);
else
    l = length(kernel);
    I = eye(l);
    alpha = 2;
    w = chebfun(@(x)(1-x^2)^(alpha/2));
    U = w.*chebfun(jac2cheb(I, alpha, alpha), 'coeffs');
    % Normalize
    List_I = (0:l-1);
    P_l = 1;
    for i = 1:alpha
        P_l = P_l.*((List_I+i)./(List_I+alpha+i));
    end
    N = 2^(2*alpha+1)*P_l./(2*List_I+2*alpha+1);
    U = U./sqrt(N);
    S = diag(kernel);
    L = U*sqrt(S);
end

% Generate a Gaussian matrix of size l x n
u = randn(l,n);
if size(L,2) > l
    L = L';
end

% Sample the Gaussian process
f = L*u;
end
