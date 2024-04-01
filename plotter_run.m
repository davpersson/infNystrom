% Run this script to plot the results
clc
clear
addpath('Other')
addpath('results')


for test_case = 1:5
    
    if test_case == 1
        
        filename = 'results/pretty_function';
        title_for_plot = 'Exact kernel';
        plot_type = 1;
        
    elseif test_case == 2
        
        filename = 'results/matern32';
        title_for_plot = '$G_{3/2}$';
        plot_type = 2;
        
    elseif test_case == 3
        
        filename = 'results/matern52';
        title_for_plot = '$G_{5/2}$';
        plot_type = 2;
        
    elseif test_case == 4
        
        filename = 'results/matern12';
        title_for_plot = '$G_{1/2}$';
        plot_type = 2;
        
    elseif test_case == 5
        
        filename = 'results/ode';
        title_for_plot = 'Green''s functions';
        plot_type = 2;
        
    end
    
    plotter(filename,plot_type,title_for_plot)
end