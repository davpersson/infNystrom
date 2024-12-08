function test(Afun,F,U_true,S_true,covariance_cell,g,rank_list,filename)

error_list = zeros(length(covariance_cell),length(rank_list));
error_optimal = zeros(1,length(rank_list));

rank_counter = 0;
for r = rank_list
    
    rank_counter = rank_counter + 1;
    
    error_optimal(rank_counter) = norm(S_true((r+1):end,(r+1):end),'fro');
    approx_cell = cell(1,length(covariance_cell));
    
    for i = 1:length(covariance_cell)
        
        %Unpack covariance cell
        C = covariance_cell{1,i};
        V = C{1,1}; V = V{1,1}; Lambda = C{1,2}; Lambda = Lambda{1,1};
        
        % Run the Nyström approximation
        [U,S] = nystrom(Afun,V,Lambda,r);
        approx_cell{1,i} = {{U},{S}};
        
        %Compute error in the Hilbert-Schmidt norm
        %inprod = U_true'*U;
        gS_true = diag(g(diag(S_true)));
        gS = diag(g(diag(S)));
        %error_list(i,rank_counter) = sqrt(norm(gS_true,'fro')^2 + norm(gS,'fro')^2-2*trace(gS_true*inprod*gS*inprod'));
        error_list(i,rank_counter) = norm(U_true*gS_true*U_true' - U*gS*U','fro');
        
    end
    
end

save(filename,'error_list','rank_list','error_optimal','approx_cell','U_true','S_true','g')

end