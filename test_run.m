% Run this script to run the tests

clc
clear

rng(0)
addpath('Other')
addpath('results')
tic

for test_case = 1:6
    
    test_case
    
    if test_case == 1
        
        f = @(x,y) 1./(1+100*(x.^2-y.^2).^2);
        box = [-1 1 -1 1];
        F = chebfun2(f);
        [U_true,S_true,~] = svd(F);
        Afun = @(X) U_true*(S_true*(U_true'*X));
        rank_list = 5:5:50;
        g = @(x) x;
        filename = 'results/pretty_function';
        
        %k = @(x,y,l) exp(-((x-y).^2)/(2*l^2))./sqrt(2*pi*l^2);
        k = @(x,y,l) exp(-((x-y).^2)/(2*l^2));
        [V1,Lambda1,~] = svd(chebfun2(@(x,y) k(x,y,1)));
        [V2,Lambda2,~] = svd(chebfun2(@(x,y) k(x,y,0.1)));
        [V3,Lambda3,~] = svd(chebfun2(@(x,y) k(x,y,0.01)));
        covariance_cell = {{{V1},{Lambda1}},{{V2},{Lambda2}},{{V3},{Lambda3}}};
        %covariance_cell = {{{U_true},{S_true}},{{U_true},{S_true}},{{U_true},{S_true}}};

        
    elseif test_case == 2
        
        f = @(x,y) (1+sqrt(3)*abs(x-y)).*exp(-sqrt(3)*abs(x-y));
        box = [-1 1 -1 1];
        F = chebfun2(f);
        [U_true,S_true,~] = svd(F);
        Afun = @(X) U_true*(S_true*(U_true'*X));
        rank_list = 10:10:100;
        g = @(x) x;
        filename = 'results/matern32';
        
        %k = @(x,y,l) exp(-((x-y).^2)/(2*l^2))./sqrt(2*pi*l^2);
        k = @(x,y,l) exp(-((x-y).^2)/(2*l^2));
        [V1,Lambda1,~] = svd(chebfun2(@(x,y) k(x,y,1)));
        [V2,Lambda2,~] = svd(chebfun2(@(x,y) k(x,y,0.1)));
        [V3,Lambda3,~] = svd(chebfun2(@(x,y) k(x,y,0.01)));
        covariance_cell = {{{V1},{Lambda1}},{{V2},{Lambda2}},{{V3},{Lambda3}}};
        
    elseif test_case == 3
        
        f = @(x,y) (1+sqrt(5)*abs(x-y) + (5/3)*(x-y).^2).*exp(-sqrt(5)*abs(x-y));
        box = [-1 1 -1 1];
        F = chebfun2(f);
        [U_true,S_true,~] = svd(F);
        Afun = @(X) U_true*(S_true*(U_true'*X));
        rank_list = 3:3:30;
        g = @(x) x;
        filename = 'results/matern52';
        
        %k = @(x,y,l) exp(-((x-y).^2)/(2*l^2))./sqrt(2*pi*l^2);
        k = @(x,y,l) exp(-((x-y).^2)/(2*l^2));
        [V1,Lambda1,~] = svd(chebfun2(@(x,y) k(x,y,1)));
        [V2,Lambda2,~] = svd(chebfun2(@(x,y) k(x,y,0.1)));
        [V3,Lambda3,~] = svd(chebfun2(@(x,y) k(x,y,0.01)));
        covariance_cell = {{{V1},{Lambda1}},{{V2},{Lambda2}},{{V3},{Lambda3}}};
        
    elseif test_case == 4
        
        f = @(x,y) exp(-abs(x-y));
        box = [-1 1 -1 1];
        F = chebfun2(f);
        [U_true,S_true,~] = svd(F);
        Afun = @(X) U_true*(S_true*(U_true'*X));
        rank_list = 10:10:100;
        g = @(x) x;
        filename = 'results/matern12';
        
        %k = @(x,y,l) exp(-((x-y).^2)/(2*l^2))./sqrt(2*pi*l^2);
        k = @(x,y,l) exp(-((x-y).^2)/(2*l^2));
        [V1,Lambda1,~] = svd(chebfun2(@(x,y) k(x,y,1)));
        [V2,Lambda2,~] = svd(chebfun2(@(x,y) k(x,y,0.1)));
        [V3,Lambda3,~] = svd(chebfun2(@(x,y) k(x,y,0.01)));
        covariance_cell = {{{V1},{Lambda1}},{{V2},{Lambda2}},{{V3},{Lambda3}}};
        
    elseif test_case == 5
        
        f = @(x,y) min(x,y);
        box = [0 1 0 1];
        F = chebfun2(f,box);
        [U_true,S_true,~] = svd(F);
        Afun = @(X) U_true*(S_true*(U_true'*X));
        rank_list = 10:10:100;
        g = @(x) x;
        filename = 'results/brownian_motion';
        
        %k = @(x,y,l) exp(-((x-y).^2)/(2*l^2))./sqrt(2*pi*l^2);
        k = @(x,y,l) exp(-((x-y).^2)/(2*l^2));
        [V1,Lambda1,~] = svd(chebfun2(@(x,y) k(2*x-1,2*y-1,1),box));
        [V2,Lambda2,~] = svd(chebfun2(@(x,y) k(2*x-1,2*y-1,0.1),box));
        [V3,Lambda3,~] = svd(chebfun2(@(x,y) k(2*x-1,2*y-1,0.01),box));
        covariance_cell = {{{V1},{Lambda1}},{{V2},{Lambda2}},{{V3},{Lambda3}}};
        
    elseif test_case == 6
        
        f = @(x,y) min(x,y) - x.*y/(2*pi);
        box = [0 2*pi 0 2*pi];
        F = chebfun2(f,box);
        [U_true,S_true,~] = svd(F);
        Afun = @(X) U_true*(S_true*(U_true'*X));
        rank_list = 10:10:100;
        g = @(x) x./(x+1);
        filename = 'results/ode';
        
        %k = @(x,y,l) exp(-((x-y).^2)/(2*l^2))./sqrt(2*pi*l^2);
        k = @(x,y,l) exp(-((x-y).^2)/(2*l^2));
        [V1,Lambda1,~] = svd(chebfun2(@(x,y) k(x/pi-1,y/pi-1,1),box));
        [V2,Lambda2,~] = svd(chebfun2(@(x,y) k(x/pi-1,y/pi-1,0.1),box));
        [V3,Lambda3,~] = svd(chebfun2(@(x,y) k(x/pi-1,y/pi-1,0.01),box));
        covariance_cell = {{{V1},{Lambda1}},{{V2},{Lambda2}},{{V3},{Lambda3}}};
        
    end
    
    %Run test
    test(Afun,F,U_true,S_true,covariance_cell,g,rank_list,filename)
    
    %Plot results
    plotter(filename)
    
end
toc