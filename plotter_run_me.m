% Run this script to create the figures. Run run_me.m to rerun all experiments 

clc
clear

rng(0)
addpath('Other')
addpath('results')

tic

plotter_pretty_function_example('results/pretty_function') % Figure 1
plotter_matern_example('results/matern12','results/matern52',0) % Figure 2
plotter_gaussian_process_example('results/gaussian_process') % Figure 3
plotter_inverse_example('results/inverse_problem') % Figure 4

toc