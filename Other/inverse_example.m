function inverse_example()
% This function creates Figure 5

% Set up experiment parameters
N = 50; % Number of nodes in discretization grid
lambda = 0.1; % Value of lambda
t_list = [0.1 0.5 1]; % List of times

if t_list~= [0.1 0.5 1]

    warning('The plotter function will assume that t_list = [0.1 0.5 1]')

end

% --- Set up experiment ---
A = uq_laplace(N,lambda,t_list); % Matrix representing the discretization of the operator
Afun = @(X) A*X; % Function to access operator through matvecs
rank_list = 1:30; % List of ranks

% Covariance operator
T = 2*eye(N) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1);
L = kron(T,eye(N)) + kron(eye(N),T);
[V,Lambda,~] = svd(L);
Lambda = diag(diag(Lambda).^(-2)); 

filename = 'results/inverse_problem'; % Name of the file storing the results
% ------

s_true = sort(svd(A),'descend'); % List of singular values in descending order
error_optimal = sum(s_true) - cumsum(s_true); % Trace errors
error_optimal = error_optimal(1:rank_list(end)); % Optimal low-rank approximation errors
error_list = zeros(1,length(rank_list));

rank_counter = 0; % Iteration counter
for r = rank_list
    
    rank_counter = rank_counter + 1 %Update counter

    [~,S] = nystrom(Afun,V,Lambda,r); % Run the Nyström approximation
    error_list(rank_counter) = sum(s_true) - trace(S); % Compute errors
    
end
%------

% --- Get a few solutions ---
solution = uq_laplace_sol(N,lambda,t_list);

% --- Save results ---
save(filename,'error_list','rank_list','error_optimal','s_true','solution','N');
% ------

% --- Plot results ---
plotter_inverse_example(filename)
% ------

end

function solution = uq_laplace_sol(N,lambda,t_list)
% This functions solves the PDE for the times indicated in t_list. N is the
% number of discretization points

h = 2/(N+1);

T = -(2*eye(N) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1)) / h^2;
L = kron(T,eye(N)) + kron(eye(N),T);

[U,Lambda] = eig(L + lambda*eye(N^2));

F = eye(N^2);

for t = t_list
    
    F = [F; U*diag(exp(t*diag(Lambda)))*U'];
    
end

initial_condition = L\randn(N^2,1);

solution_temp = F*initial_condition;
solution = zeros(N,N, length(t_list)+1);

for k = 1:(length(t_list) + 1)

    solution(:,:,k) = reshape(solution_temp(((k-1)*N^2+1):k*N^2,:),[N N]);

end

end

function A = uq_laplace(N,lambda,t_list)
% Function creating a discretization of the operator

h = 2/(N+1);

T = -(2*eye(N) - diag(ones(N-1,1),1) - diag(ones(N-1,1),-1)) / h^2;
L = kron(T,eye(N)) + kron(eye(N),T);

[U,Lambda] = eig(L + lambda*eye(N^2));

F = [];

for t = t_list
    
    F = [F; U*diag(exp(t*diag(Lambda)))*U'];
    
end

A = F/L;
A = A'*A;
end