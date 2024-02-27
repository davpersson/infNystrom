% Run this script to plot the results
clc
clear


for test_case = 1:6
    
    if test_case == 1
        
        filename = 'results/pretty_function';
        
    elseif test_case == 2
        
        filename = 'results/matern32';
        
    elseif test_case == 3
        
        filename = 'results/matern52';
        
    elseif test_case == 4
        
        filename = 'results/matern12';
        
    elseif test_case == 5
        
        filename = 'results/brownian_motion';
        
    elseif test_case == 6
        
        filename = 'results/ode';
        
    end
    
    plotter(filename)
end