G = chebfun2(@(x,y)min(x,y)-x.*y/(2*pi),[0,1,0,1],[500,500]);
[U,S,~] = svd(G);
l = 0.01;
gaussian_01 = chebfun2(@(x,y) exp(-(x-y).^2/(2*l^2)),[0,1,0,1]);

K = U'*gaussian_01*U;
k = 50;
K11 = K(1:k,1:k);
K21 = K(k+1:end,1:k);
K22 = K(k+1:end,k+1:end);
S1 = S(1:k,1:k);
S2 = S(k+1:end,k+1:end);
K221 = K22 - (K21 / K11)*K21';
beta_k = norm(S2.^0.5*K221*S2.^0.5,'fro')*norm(inv(K11),2)/norm(S2,'fro')
delta_k = norm(S2.^0.5*(K21/K11)*K21'*S2.^0.5,'fro')*norm(inv(K11),2)/norm(S2,'fro')