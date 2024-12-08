function Fn = randomsvd(F, varargin)
% RANDOMSVD 
%
% RANDOMSVD(f, K, n) returns a rank n approximation to a self-adjoint linear chebop
% operator using n samples from a Gaussian process with covariance kernel K
%
% RANDOMSVD(f, f_adjoint K, n) returns a rank n approximation to a linear chebop
% operator F with adjoint f_adjoint using n samples from a Gaussian process with covariance kernel K

% Parse the input
n = varargin{end};
kernel = varargin{end-1};
if nargin > 3
    F_adjoint = varargin{1};
else
    F_adjoint = F;
end

% Sample the Gaussian process
Omega = gaussianprocess(kernel, n);
A = F(Omega);

% Compute rank n chebfun2 approximation to the kernel
Q = orth(A);
% Use symmetry of the integral kernel to compute QQ'*F = Q*(F(Q))'
Fn = F_adjoint(Q)*Q'; % Need to transpose to recover the kernel
end