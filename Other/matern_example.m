function matern_example()
% This function runs the numerical experiments for Figure 2 and Figure 3.

% --- Set up experiment ---
f1 = @(x,y) exp(-abs(x-y)); % Kernel defining the integral operator (Matern-1/2 kernel)
F1 = chebfun2(f1,[500,500]); % Create a chebfun2 object
[U_true1,S_true1,~] = svd(F1); % Compute the SVD of the integral operator
Afun1 = @(X) U_true1*(S_true1*(U_true1'*X)); % Function to access operator through matvecs
filename1 = 'results/matern12'; % Name of the file storing the results

f2 = @(x,y) (1+sqrt(5)*abs(x-y) + (5/3)*(x-y).^2).*exp(-sqrt(5)*abs(x-y)); % Kernel defining the integral operator (Matern-5/2 kernel)
F2 = chebfun2(f2,[500 500]); % Create a chebfun2 object
[U_true2,S_true2,~] = svd(F2); % Compute the SVD of the integral operator
Afun2 = @(X) U_true2*(S_true2*(U_true2'*X)); % Function to access operator through matvecs
filename2 = 'results/matern52'; % Name of the file storing the results

rank_list = 1:1:100; % List of all ranks

% Covariance operators 
k = @(x,y,l) exp(-((x-y).^2)/(2*l^2)); % Squared-exponential kernel with length parameter l
[V1,Lambda1,~] = svd(chebfun2(@(x,y) k(x,y,1)));
[V2,Lambda2,~] = svd(chebfun2(@(x,y) (1+sqrt(3)*abs(x-y)).*exp(-sqrt(3)*abs(x-y)),[500 500])); % Matern-3/2 kernel
[V3,Lambda3,~] = svd(chebfun2(@(x,y) k(x,y,0.01)));
covariance_cell = {{{V1},{Lambda1}},{{V2},{Lambda2}},{{V3},{Lambda3}}}; % A cell containing the Mercer expansions of the covariance kernels for the random fields
%------

% --- Run tests and plot the results ---
matern_test(Afun1,U_true1,S_true1,covariance_cell,rank_list,filename1,F1);
matern_test(Afun2,U_true2,S_true2,covariance_cell,rank_list,filename2,F2);
plotter_matern_example(filename1,filename2,1)
%------
end


function matern_test(Afun,U_true,S_true,covariance_cell,rank_list,filename,F)
% Function  to run the numerical experiments in Figure 1

error_list = zeros(length(covariance_cell),length(rank_list)); %Allocate space for errors

% Compute optimal errors
error_optimal = sum(diag(S_true)) - cumsum(diag(S_true));
error_optimal = error_optimal(1:rank_list(end));

rank_counter = 0; % Counter
for r = rank_list
    
    rank_counter = rank_counter + 1 % Update counter
        
    for i = 1:length(covariance_cell)
        
        %Unpack covariance cell
        C = covariance_cell{1,i};
        V = C{1,1}; V = V{1,1}; Lambda = C{1,2}; Lambda = Lambda{1,1};
        
        % Run the Nyström approximation
        [U,S] = nystrom(Afun,V,Lambda,r);

        %Compute error
        error_list(i,rank_counter) = trace(S_true) - trace(S);
        
    end

end

% Save results
save(filename,'error_list','rank_list','error_optimal','U_true','S_true','covariance_cell','F')

end