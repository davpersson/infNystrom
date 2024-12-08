function pretty_function_example()
% This function runs the numerical experiments for Figure 1.

% --- Set up experiment ---
f = @(x,y) 1./(1+100*(x.^2-y.^2).^2); % Kernel defining the integral operator
F = chebfun2(f); % Create a chebfun2 object
[U_true,S_true,~] = svd(F); % Compute the SVD of the integral operator
Afun = @(X) U_true*(S_true*(U_true'*X)); % Function to access operator through matvecs
rank_list = 1:1:100; % List of all ranks
filename = 'results/pretty_function'; % Name of the file storing the results

% Covariance operators
k = @(x,y,l) exp(-((x-y).^2)/(2*l^2)); % Squared-exponential kernel with length parameter l
[V1,Lambda1,~] = svd(chebfun2(@(x,y) k(x,y,1)));
[V2,Lambda2,~] = svd(chebfun2(@(x,y) k(x,y,0.1)));
[V3,Lambda3,~] = svd(chebfun2(@(x,y) k(x,y,0.01)));
covariance_cell = {{{V1},{Lambda1}},{{V2},{Lambda2}},{{V3},{Lambda3}}}; % A cell containing the Mercer expansions of the covariance kernels for the random fields
%------

% --- Run tests and plot the results
pretty_test(Afun,U_true,S_true,covariance_cell,rank_list,filename); % Run test
plotter_pretty_function_example(filename) % Plot results
%------
end

function pretty_test(Afun,U_true,S_true,covariance_cell,rank_list,filename)
% Function  to run the numerical experiments in Figure 1

error_list = zeros(length(covariance_cell),length(rank_list)); %Allocate space for errors

% Allocate space for optimal errors
error_optimal = zeros(1,length(rank_list));

rank_counter = 0; % Counter
for r = rank_list
    
    rank_counter = rank_counter + 1 % Update counter
    error_optimal(rank_counter) = norm(S_true(r+1:end,r+1:end),'fro'); % Optimal error
    
    approx_cell = cell(1,length(covariance_cell)); % Cell containing the Nystrom approximations so that we can plot them
    
    for i = 1:length(covariance_cell)
        
        %Unpack covariance cell
        C = covariance_cell{1,i};
        V = C{1,1}; V = V{1,1}; Lambda = C{1,2}; Lambda = Lambda{1,1};
        
        % Run the Nyström approximation
        [U,S] = nystrom(Afun,V,Lambda,r);
        approx_cell{1,i} = {{U},{S}}; % Save Nystrom approximation

        %Compute error
        error_list(i,rank_counter) = norm(U_true*S_true*U_true'-U*S*U','fro');
        
    end

end

% Save results
save(filename,'error_list','rank_list','error_optimal','U_true','S_true','covariance_cell','approx_cell')

end