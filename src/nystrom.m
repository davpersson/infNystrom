function Fn = nystrom(F, K, n)
% NYSTROM 
%
% NYSTROM(F, K, n) returns a rank n approximation to a self-adjoint linear chebop
% operator using n samples from a Gaussian process with covariance kernel K

% Use Algorithm 2.1 in https://arxiv.org/pdf/2009.11392.pdf for stability

% Oversampling parameter
p = floor(0.5*n);

% Sample the Gaussian process
X = gaussianprocess(K, n);
Y = gaussianprocess(K, n+p);


% Sample A at X and Y
Ax = F*X;
Ay = F*Y;
Axy = Y'*Ax;

% Compute rank n chebfun2 approximation to the kernel
[Q, R] = qr(Axy);

% Truncate R
R = R(1:n,1:n);
Q = Q(:,1:n);

% Compute approximation
Fn = (R'\Ax')'*(Q'*Ay');
end